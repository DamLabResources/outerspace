"""Tests for the pipeline command and Snakefile location logic"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock, mock_open
from outerspace.cli.main import Cli
from outerspace.cli.commands.pipeline import PipelineCommand


@pytest.fixture
def mock_pipeline_command():
    """Fixture to create a PipelineCommand instance with mock args"""
    cmd = PipelineCommand()
    cmd.args = Mock()
    cmd.args.config_file = 'test_config.toml'
    cmd.args.snakemake_config = 'test_snakemake.yaml'
    cmd.args.snakemake_args = None
    return cmd


@pytest.fixture
def temp_workflow_dir(tmp_path):
    """Fixture to create a temporary workflow directory structure"""
    # Create workflow directory structure
    workflow_dir = tmp_path / 'workflow'
    workflow_dir.mkdir()
    
    snakefile = workflow_dir / 'Snakefile'
    snakefile.write_text('# Test Snakefile\nrule all:\n    input: "test.txt"')
    
    # Create wrappers
    for wrapper in ['findseq', 'collapse', 'count', 'merge', 'stats']:
        wrapper_dir = workflow_dir / 'wrappers' / wrapper
        wrapper_dir.mkdir(parents=True)
        (wrapper_dir / 'wrapper.py').write_text('# Test wrapper')
        (wrapper_dir / 'environment.yaml').write_text('name: test')
        (wrapper_dir / 'README.md').write_text('# Test')
    
    return workflow_dir


def test_pipeline_initialization():
    """Test that pipeline command initializes correctly"""
    args = [
        "pipeline",
        "test_config.toml",
        "test_snakemake.yaml",
    ]
    cli = Cli(args)
    assert cli.args.config_file == "test_config.toml"
    assert cli.args.snakemake_config == "test_snakemake.yaml"
    assert cli.args.snakemake_args is None


def test_pipeline_with_snakemake_args():
    """Test that pipeline command accepts additional Snakemake arguments"""
    args = [
        "pipeline",
        "test_config.toml",
        "test_snakemake.yaml",
        "--snakemake-args",
        "--dry-run --cores 4",
    ]
    cli = Cli(args)
    assert cli.args.snakemake_args == "--dry-run --cores 4"


def test_get_repo_snakefile_exists(tmp_path):
    """Test locating Snakefile in repository"""
    # Create a mock command instance
    with patch('outerspace.cli.commands.pipeline.__file__', 
               str(tmp_path / 'outerspace' / 'cli' / 'commands' / 'pipeline.py')):
        # Create workflow directory structure
        workflow_dir = tmp_path / 'workflow'
        workflow_dir.mkdir()
        snakefile = workflow_dir / 'Snakefile'
        snakefile.write_text('# Test Snakefile')
        
        # Create command instance
        cmd = PipelineCommand()
        
        # Test repo snakefile location
        result = cmd._get_repo_snakefile()
        assert result is not None
        assert result.exists()
        assert result.name == 'Snakefile'


def test_get_repo_snakefile_not_exists(tmp_path):
    """Test when repository Snakefile doesn't exist"""
    with patch('outerspace.cli.commands.pipeline.__file__',
               str(tmp_path / 'outerspace' / 'cli' / 'commands' / 'pipeline.py')):
        cmd = PipelineCommand()
        result = cmd._get_repo_snakefile()
        assert result is None


def test_locate_snakefile_from_repo(temp_workflow_dir):
    """Test that _locate_snakefile finds repository Snakefile"""
    snakefile = temp_workflow_dir / 'Snakefile'
    
    with patch('outerspace.cli.commands.pipeline.__file__',
               str(temp_workflow_dir.parent / 'outerspace' / 'cli' / 'commands' / 'pipeline.py')):
        cmd = PipelineCommand()
        
        # Mock _get_packaged_snakefile to return None
        with patch.object(cmd, '_get_packaged_snakefile', return_value=None):
            result = cmd._locate_snakefile()
            assert result is not None
            assert result.exists()
            assert result.name == 'Snakefile'


def test_locate_snakefile_not_found():
    """Test that _locate_snakefile raises error when Snakefile not found"""
    cmd = PipelineCommand()
    
    # Mock both methods to return None
    with patch.object(cmd, '_get_packaged_snakefile', return_value=None):
        with patch.object(cmd, '_get_repo_snakefile', return_value=None):
            with pytest.raises(ValueError, match="Unable to locate Snakefile"):
                cmd._locate_snakefile()


def test_locate_snakefile_prefers_packaged(tmp_path):
    """Test that _locate_snakefile prefers packaged Snakefile over repo"""
    packaged_path = tmp_path / 'packaged' / 'Snakefile'
    packaged_path.parent.mkdir(parents=True)
    packaged_path.write_text('# Packaged Snakefile')
    
    repo_path = tmp_path / 'repo' / 'Snakefile'
    repo_path.parent.mkdir(parents=True)
    repo_path.write_text('# Repo Snakefile')
    
    cmd = PipelineCommand()
    
    # Mock to return packaged path
    with patch.object(cmd, '_get_packaged_snakefile', return_value=packaged_path):
        with patch.object(cmd, '_get_repo_snakefile', return_value=repo_path):
            result = cmd._locate_snakefile()
            assert result == packaged_path


def test_build_snakemake_argv_basic(mock_pipeline_command, tmp_path):
    """Test building basic Snakemake argv"""
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    argv = mock_pipeline_command._build_snakemake_argv(snakefile, [])
    
    assert argv[0] == 'snakemake'
    assert argv[1] == '-s'
    assert argv[2] == str(snakefile)
    assert '--configfile' in argv
    assert 'test_snakemake.yaml' in argv
    assert '--config' in argv
    assert 'toml=test_config.toml' in argv


def test_build_snakemake_argv_with_user_args(mock_pipeline_command, tmp_path):
    """Test building Snakemake argv with user arguments"""
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    user_args = ['--dry-run', '--cores', '4']
    argv = mock_pipeline_command._build_snakemake_argv(snakefile, user_args)
    
    assert '--dry-run' in argv
    assert '--cores' in argv
    assert '4' in argv


def test_build_snakemake_argv_respects_user_configfile(mock_pipeline_command, tmp_path):
    """Test that user-provided --configfile is not overridden"""
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    user_args = ['--configfile', 'custom.yaml']
    argv = mock_pipeline_command._build_snakemake_argv(snakefile, user_args)
    
    # Should only have one --configfile (the user's)
    configfile_count = argv.count('--configfile')
    assert configfile_count == 1
    configfile_idx = argv.index('--configfile')
    assert argv[configfile_idx + 1] == 'custom.yaml'


def test_build_snakemake_argv_respects_user_config(mock_pipeline_command, tmp_path):
    """Test that user-provided --config is not overridden"""
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    user_args = ['--config', 'custom_key=custom_value']
    argv = mock_pipeline_command._build_snakemake_argv(snakefile, user_args)
    
    # Should only have one --config (the user's)
    config_count = argv.count('--config')
    assert config_count == 1


def test_parse_snakemake_args_empty():
    """Test parsing empty Snakemake args"""
    cmd = PipelineCommand()
    result = cmd._parse_snakemake_args(None)
    assert result == []
    
    result = cmd._parse_snakemake_args("")
    assert result == []


def test_parse_snakemake_args_simple():
    """Test parsing simple Snakemake args"""
    cmd = PipelineCommand()
    result = cmd._parse_snakemake_args("--dry-run --cores 4")
    assert result == ['--dry-run', '--cores', '4']


def test_parse_snakemake_args_with_quotes():
    """Test parsing Snakemake args with quoted values"""
    cmd = PipelineCommand()
    result = cmd._parse_snakemake_args('--config "key=value with spaces"')
    assert result == ['--config', 'key=value with spaces']


def test_execute_snakemake_success():
    """Test successful Snakemake execution"""
    cmd = PipelineCommand()
    
    with patch('outerspace.cli.commands.pipeline.snakemake.main') as mock_snakemake:
        # Mock successful execution (no exception)
        mock_snakemake.return_value = None
        
        # Should not raise
        cmd._execute_snakemake(['snakemake', '--help'])
        mock_snakemake.assert_called_once_with(['--help'])


def test_execute_snakemake_failure():
    """Test failed Snakemake execution"""
    cmd = PipelineCommand()
    
    with patch('outerspace.cli.commands.pipeline.snakemake.main') as mock_snakemake:
        # Mock failed execution
        mock_snakemake.side_effect = SystemExit(1)
        
        with pytest.raises(SystemExit):
            cmd._execute_snakemake(['snakemake', '--help'])


def test_execute_snakemake_success_exit():
    """Test Snakemake execution with SystemExit(0)"""
    cmd = PipelineCommand()
    
    with patch('outerspace.cli.commands.pipeline.snakemake.main') as mock_snakemake:
        # Mock successful exit
        mock_snakemake.side_effect = SystemExit(0)
        
        # Should not raise
        cmd._execute_snakemake(['snakemake', '--help'])


@pytest.mark.integration
def test_get_packaged_snakefile_python39_style(tmp_path):
    """Test getting packaged Snakefile using Python 3.9+ importlib.resources"""
    cmd = PipelineCommand()
    
    # Mock the resource_files approach
    mock_resource = Mock()
    mock_resource.is_file.return_value = True
    mock_snakefile_path = tmp_path / 'Snakefile'
    mock_snakefile_path.write_text('# Test')
    
    with patch('outerspace.cli.commands.pipeline._HAS_FILES_API', True):
        with patch('outerspace.cli.commands.pipeline.resource_files') as mock_res_files:
            mock_res_files.return_value.joinpath.return_value = mock_resource
            with patch('pathlib.Path.__new__') as mock_path:
                mock_path.return_value = mock_snakefile_path
                # This would test the Python 3.9+ code path
                # Note: Full test requires actual importlib.resources behavior


@pytest.mark.integration  
def test_get_packaged_snakefile_python38_style(tmp_path):
    """Test getting packaged Snakefile using Python 3.8 fallback"""
    # Create a mock outerspace package structure
    pkg_dir = tmp_path / 'site-packages' / 'outerspace'
    pkg_dir.mkdir(parents=True)
    (pkg_dir / '__init__.py').write_text('')
    
    workflow_dir = pkg_dir / 'workflow'
    workflow_dir.mkdir()
    snakefile = workflow_dir / 'Snakefile'
    snakefile.write_text('# Test Snakefile')
    
    cmd = PipelineCommand()
    
    # Mock the Python 3.8 fallback
    with patch('outerspace.cli.commands.pipeline._HAS_FILES_API', False):
        mock_outerspace = Mock()
        mock_outerspace.__file__ = str(pkg_dir / '__init__.py')
        
        with patch('importlib.import_module', return_value=mock_outerspace):
            result = cmd._get_packaged_snakefile()
            # In real scenario, this should find the Snakefile


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.

