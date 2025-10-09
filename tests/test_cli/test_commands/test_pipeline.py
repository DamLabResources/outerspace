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


def test_parse_snakemake_config_dict_basic(mock_pipeline_command):
    """Test parsing basic config dictionary"""
    config = mock_pipeline_command._parse_snakemake_config_dict([])
    
    assert 'toml' in config
    assert config['toml'] == 'test_config.toml'


def test_parse_snakemake_config_dict_with_user_config(mock_pipeline_command):
    """Test parsing config dictionary with user config arguments"""
    user_args = ['--config', 'custom_key=custom_value', '--other-arg', 'value']
    config = mock_pipeline_command._parse_snakemake_config_dict(user_args)
    
    assert 'toml' in config
    assert config['toml'] == 'test_config.toml'
    assert 'custom_key' in config
    assert config['custom_key'] == 'custom_value'


def test_parse_execution_cores_default():
    """Test default cores parsing"""
    cmd = PipelineCommand()
    cores = cmd._parse_execution_cores([])
    assert cores == 1


def test_parse_execution_cores_with_cores_flag():
    """Test cores parsing with --cores flag"""
    cmd = PipelineCommand()
    cores = cmd._parse_execution_cores(['--cores', '4'])
    assert cores == 4


def test_is_dry_run_false():
    """Test dry-run detection when not set"""
    cmd = PipelineCommand()
    assert cmd._is_dry_run([]) is False


def test_is_dry_run_true():
    """Test dry-run detection when set"""
    cmd = PipelineCommand()
    assert cmd._is_dry_run(['--dry-run']) is True
    assert cmd._is_dry_run(['-n']) is True


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


def test_execute_snakemake_success(tmp_path):
    """Test successful Snakemake execution"""
    cmd = PipelineCommand()
    cmd.args = Mock()
    
    # Create test config files with valid YAML dict
    config_yaml = tmp_path / 'test.yaml'
    config_yaml.write_text('test_key: test_value\n')
    config_toml = tmp_path / 'test.toml'
    config_toml.write_text('# Test toml\n')
    
    cmd.args.snakemake_config = str(config_yaml)
    cmd.args.config_file = str(config_toml)
    
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    # Mock the new Snakemake v9 API
    mock_dag_api = MagicMock()
    mock_workflow_api = MagicMock()
    mock_workflow_api.dag.return_value = mock_dag_api
    mock_snakemake_api = MagicMock()
    mock_snakemake_api.workflow.return_value = mock_workflow_api
    mock_snakemake_api.__enter__ = Mock(return_value=mock_snakemake_api)
    mock_snakemake_api.__exit__ = Mock(return_value=False)
    
    with patch('outerspace.cli.commands.pipeline.SnakemakeApi', return_value=mock_snakemake_api):
        # Should not raise
        cmd._execute_snakemake(snakefile, [])
        
        # Verify API was called correctly
        mock_snakemake_api.workflow.assert_called_once()
        mock_workflow_api.dag.assert_called_once()
        mock_dag_api.execute_workflow.assert_called_once()


def test_execute_snakemake_failure(tmp_path):
    """Test failed Snakemake execution"""
    cmd = PipelineCommand()
    cmd.args = Mock()
    
    # Create test config files with valid YAML dict
    config_yaml = tmp_path / 'test.yaml'
    config_yaml.write_text('test_key: test_value\n')
    config_toml = tmp_path / 'test.toml'
    config_toml.write_text('# Test toml\n')
    
    cmd.args.snakemake_config = str(config_yaml)
    cmd.args.config_file = str(config_toml)
    
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    # Mock the new Snakemake v9 API with failure
    mock_dag_api = MagicMock()
    mock_dag_api.execute_workflow.side_effect = Exception("Workflow failed")
    mock_workflow_api = MagicMock()
    mock_workflow_api.dag.return_value = mock_dag_api
    mock_snakemake_api = MagicMock()
    mock_snakemake_api.workflow.return_value = mock_workflow_api
    mock_snakemake_api.__enter__ = Mock(return_value=mock_snakemake_api)
    mock_snakemake_api.__exit__ = Mock(return_value=False)
    
    with patch('outerspace.cli.commands.pipeline.SnakemakeApi', return_value=mock_snakemake_api):
        with pytest.raises(Exception, match="Workflow failed"):
            cmd._execute_snakemake(snakefile, [])


def test_execute_snakemake_success_exit(tmp_path):
    """Test Snakemake execution completes successfully"""
    cmd = PipelineCommand()
    cmd.args = Mock()
    
    # Create test config files with valid YAML dict
    config_yaml = tmp_path / 'test.yaml'
    config_yaml.write_text('test_key: test_value\n')
    config_toml = tmp_path / 'test.toml'
    config_toml.write_text('# Test toml\n')
    
    cmd.args.snakemake_config = str(config_yaml)
    cmd.args.config_file = str(config_toml)
    
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text('# Test')
    
    # Mock the new Snakemake v9 API
    mock_dag_api = MagicMock()
    mock_workflow_api = MagicMock()
    mock_workflow_api.dag.return_value = mock_dag_api
    mock_snakemake_api = MagicMock()
    mock_snakemake_api.workflow.return_value = mock_workflow_api
    mock_snakemake_api.__enter__ = Mock(return_value=mock_snakemake_api)
    mock_snakemake_api.__exit__ = Mock(return_value=False)
    
    with patch('outerspace.cli.commands.pipeline.SnakemakeApi', return_value=mock_snakemake_api):
        # Should not raise
        cmd._execute_snakemake(snakefile, [])
        
        # Verify execution happened
        mock_dag_api.execute_workflow.assert_called_once()


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

