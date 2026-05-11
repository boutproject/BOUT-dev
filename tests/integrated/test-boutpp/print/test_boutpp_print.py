import subprocess


def test_boutpp_print(run_isolated):

    if run_isolated():
        return

    msg = "Can we print to the log from python? 🎉 __does_it_still__work {} {:s}"
    cmd = ["python3", "test.py", msg]

    print(f"+ python3 test.py {msg} > out.log")
    with open("out.log", "w") as f:
        subprocess.run(cmd, stdout=f, check=True)

    with open("out.log", "r") as f:
        assert msg in f.read(), f"Error: '{msg}' not found in out.log"

    with open("test/BOUT.log.0", "r") as f:
        assert msg in f.read(), f"Error: '{msg}' not found in test/BOUT.log.0"
