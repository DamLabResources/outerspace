"""Main CLI entry point for OUTERSPACE.

This module provides the main command-line interface for OUTERSPACE, handling
argument parsing, command routing, and execution. It supports multiple subcommands
for different analysis workflows.
"""

import logging
from argparse import ArgumentParser, Namespace
import sys
from typing import Optional, Union

from outerspace.cli.commands.findseq import FindSeqCommand
from outerspace.cli.commands.collapse import CollapseCommand
from outerspace.cli.commands.count import CountCommand
from outerspace.cli.commands.merge import MergeCommand
from outerspace.cli.commands.stats import StatsCommand
from outerspace.cli.commands.subsample import SubsampleCommand
from outerspace.cli.commands.visualize import VisualizeCommand
from outerspace.cli.commands.pipeline import PipelineCommand

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class Cli:
    """Command-line interface for OUTERSPACE.

    This class handles argument parsing, command routing, and execution for
    the OUTERSPACE command-line tool. It supports multiple subcommands for
    different analysis workflows including sequence finding, UMI collapsing,
    counting, merging, statistics, visualization, and pipeline execution.
    """

    def __init__(self, args: Optional[Union[list, Namespace]] = None) -> None:
        """Initialize the CLI with arguments.

        Parameters
        ----------
        args : Optional[Union[list, Namespace]], default=None
            Command-line arguments. Can be a list of strings, a Namespace object,
            or None (uses sys.argv)
        """
        self.parser = self._init_parser()

        if isinstance(args, Namespace):
            self.args = args
        else:
            self.args = self.parser.parse_args(args)

        self.command = self._init_command()

        logger.debug(
            f"CLI initialized with command: {getattr(self.args, 'command', None)}"
        )

    @staticmethod
    def _init_parser() -> ArgumentParser:
        """Initialize the main argument parser with subcommands.

        Returns
        -------
        ArgumentParser
            Configured argument parser with all available subcommands

        Notes
        -----
        This method creates the main parser and registers all available commands.
        Each command is responsible for adding its own subparser and arguments.
        """
        parser = ArgumentParser(
            prog="outerspace",
            description="OUTERSPACE - Optimized Utilities for Tracking Enrichment in Screens",
            epilog="Created by ThreeBlindMice - See how they run code",
        )
        subparsers = parser.add_subparsers(dest="command", help="Available commands")

        # Register each command
        command_classes = [
            FindSeqCommand,
            CollapseCommand,
            CountCommand,
            MergeCommand,
            StatsCommand,
            SubsampleCommand,
            VisualizeCommand,
            PipelineCommand,
        ]

        for cmd_cls in command_classes:
            cmd = cmd_cls()
            cmd._init_parser(subparsers)
            logger.debug(f"Registered command: {cmd_cls.__name__}")

        return parser

    def _init_command(self) -> Optional[object]:
        """Initialize the selected command.

        Returns
        -------
        Optional[object]
            Initialized command object, or None if no command selected

        Notes
        -----
        This method maps command names to their corresponding command classes
        and initializes the selected command with the parsed arguments.
        """
        if not self.args.command:
            logger.debug("No command selected")
            return None

        command_map = {
            "findseq": FindSeqCommand,
            "collapse": CollapseCommand,
            "count": CountCommand,
            "merge": MergeCommand,
            "visualize": VisualizeCommand,
            "pipeline": PipelineCommand,
            "stats": StatsCommand,
            "subsample": SubsampleCommand,
        }

        if self.args.command not in command_map:
            logger.error(f"Unknown command: {self.args.command}")
            return None

        try:
            command = command_map[self.args.command](self.args)
            logger.debug(f"Initialized command: {self.args.command}")
            return command
        except Exception as e:
            logger.error(f"Failed to initialize command {self.args.command}: {e}")
            raise

    def run(self) -> None:
        """Run the selected command.

        If no command is selected, prints help information. Otherwise,
        executes the selected command with the parsed arguments.

        Notes
        -----
        This method handles the main execution flow of the CLI application.
        """
        if not self.command:
            logger.info("No command selected, displaying help")
            self.parser.print_help()
            return

        try:
            logger.info(f"Executing command: {self.args.command}")
            self.command.run()
            logger.info(f"Command {self.args.command} completed successfully")
        except Exception as e:
            logger.error(f"Command {self.args.command} failed: {e}")
            raise


def main(args: Optional[Union[list, Namespace]] = None) -> None:
    """Main entry point for the OUTERSPACE CLI.

    This function creates a CLI instance and runs the selected command.
    It serves as the primary entry point for the command-line application.

    Parameters
    ----------
    args : Optional[Union[list, Namespace]], default=None
        Command-line arguments. Can be a list of strings, a Namespace object,
        or None (uses sys.argv)

    Examples
    --------
    Basic usage:
        main()

    With custom arguments:
        main(['findseq', '--input', 'data.fastq'])

    From another script:
        main(Namespace(command='count', input='data.csv'))
    """
    try:
        cli = Cli(args)
        cli.run()
    except KeyboardInterrupt:
        logger.info("Command interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"CLI execution failed: {e}")
        sys.exit(1)


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
