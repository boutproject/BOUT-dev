#!/usr/bin/env python3

from boututils.run_wrapper import launch_safe
from pathlib import Path
import re
import shutil
import unittest
import platform


OUTPUT_FILE = Path("stdout.log") if platform.system() == "FreeBSD" else "stderr.log"

class TestCommandLineArgs(unittest.TestCase):
    command = "./command-args >stdout.log 2>stderr.log"

    def makeDirAndCopyInput(self, path):
        Path.mkdir(path)
        shutil.copy("BOUT.inp", path)

    def setUp(self):
        try:
            Path.unlink(OUTPUT_FILE)
        except OSError:
            pass
        shutil.rmtree("./data", ignore_errors=True)
        shutil.rmtree("./test", ignore_errors=True)
        shutil.rmtree("./test_copy", ignore_errors=True)

    # Teardown same as setup in case something went wrong with last run
    tearDown = setUp

    def testNoArgumentsNoDirectory(self):
        with self.assertRaises(RuntimeError):
            launch_safe(self.command, pipe=True, nproc=1, mthread=1)
        with open(OUTPUT_FILE) as f:
            contents = f.read()
            self.assertIn(
                '"data" does not exist',
                contents,
                msg="FAIL: Error message not printed when missing input directory",
            )

    def testHelpArgument(self):
        self.makeDirAndCopyInput("data")
        _, out = launch_safe(self.command + " --help", pipe=True, nproc=1, mthread=1)
        # Not a great test!
        self.assertNotIn(
            "Writing options", out, msg="FAIL: Attempting to write options"
        )

    def testNoArgumentsDefaultDirectory(self):
        self.makeDirAndCopyInput("data")
        launch_safe(self.command, pipe=True, nproc=1, mthread=1)
        self.assertTrue(
            Path.exists(Path("data/BOUT.settings")),
            msg="FAIL: No BOUT.settings file in data directory",
        )
        self.assertTrue(
            Path.exists(Path("data/BOUT.log.0")),
            msg="FAIL: No BOUT.log.0 file in data directory",
        )

    def testShortLogArgument(self):
        self.makeDirAndCopyInput("data")
        launch_safe(self.command + " -l different.log", pipe=True, nproc=1, mthread=1)
        self.assertFalse(
            Path.exists(Path("data/BOUT.log.0")),
            msg="FAIL: BOUT.log.0 file in data directory",
        )
        self.assertTrue(
            Path.exists(Path("data/different.log.0")),
            msg="FAIL: no different.log.0 file in data directory",
        )

    def testLongLogArgument(self):
        self.makeDirAndCopyInput("data")
        launch_safe(self.command + " --log log", pipe=True, nproc=1, mthread=1)
        self.assertFalse(
            Path.exists(Path("data/BOUT.log.0")),
            msg="FAIL: BOUT.log.0 file in data directory",
        )
        self.assertTrue(
            Path.exists(Path("data/log.0")), msg="FAIL: no log.0 file in data directory"
        )

    def testDirectoryArgument(self):
        self.makeDirAndCopyInput("test")
        launch_safe(self.command + " -d test", pipe=True, nproc=1, mthread=1)
        self.assertTrue(
            Path.exists(Path("test/BOUT.settings")),
            msg="FAIL: No BOUT.settings file in test directory",
        )
        self.assertTrue(
            Path.exists(Path("test/BOUT.log.0")),
            msg="FAIL: No BOUT.log.0 file in data directory",
        )

    def testDirectoryArgumentNonExistentDirectory(self):
        with self.assertRaises(RuntimeError):
            launch_safe(
                self.command + " -d non_existent", pipe=True, nproc=1, mthread=1
            )
        with open(OUTPUT_FILE) as f:
            contents = f.read()
            self.assertIn(
                '"non_existent" does not exist',
                contents,
                msg="FAIL: Error message not printed when missing input directory",
            )

    def testDirectoryArgumentNonDirectory(self):
        with self.assertRaises(RuntimeError):
            launch_safe(self.command + " -d runtest", pipe=True, nproc=1, mthread=1)
        with open(OUTPUT_FILE) as f:
            contents = f.read()
            self.assertIn(
                'DataDir "runtest" does not exist or is not accessible',
                contents,
                msg="FAIL: Error message not printed when missing input directory",
            )

    def testDirectoryArgumentOldSettingsFile(self):
        self.makeDirAndCopyInput("test")
        launch_safe(self.command + " -d test", pipe=True, nproc=1, mthread=1)
        shutil.copytree("test", "test_copy")
        launch_safe(
            self.command + " -d test_copy -f BOUT.settings -o testsettings",
            pipe=True,
            nproc=1,
            mthread=1,
        )

        with open("test_copy/testsettings") as f:
            contents = f.readlines()

            datadir = ""
            for line in contents:
                if line.startswith("datadir"):
                    # Assuming that this line looks like "datadir = <something> # comment"
                    datadir = line.split()[2]
                    break

            self.assertIn(
                "test_copy",
                datadir,
                msg="FAIL: datadir from command line clobbered by BOUT.settings",
            )

    def testShortOptionsAreUsed(self):
        self.makeDirAndCopyInput("test")
        launch_safe(self.command + " -d test", pipe=True, nproc=1, mthread=1)
        shutil.copytree("test", "test_copy")
        launch_safe(
            self.command + " -d test_copy -f BOUT.settings -o testsettings",
            pipe=True,
            nproc=1,
            mthread=1,
        )

        with open("test_copy/testsettings") as f:
            contents = f.read()
            matches = re.findall("not used.*Command line", contents)
            self.assertEqual(matches, [], msg="FAIL: command line options not used")

    def testCommandLineOptionsArePrinted(self):
        self.makeDirAndCopyInput("test")
        launch_safe(self.command + " -d test", pipe=True, nproc=1, mthread=1)
        shutil.copytree("test", "test_copy")
        extra_options = [
            "-d",
            "test_copy",
            "-f",
            "BOUT.settings",
            "-o",
            "testsettings",
            "-l",
            "bout.log",
            "--foo_flag",
            "some:option=value",
        ]
        launch_safe(
            self.command + " " + " ".join(extra_options), pipe=True, nproc=1, mthread=1
        )

        with open("stdout.log") as f:
            contents = f.read()

            command_line_options = None
            for line in contents.splitlines():
                if line.lstrip().startswith("Command line"):
                    command_line_options = line

            for option in extra_options:
                self.assertIn(
                    option,
                    command_line_options,
                    msg="FAIL: option {} not printed out".format(option),
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
