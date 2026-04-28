import pytest

@pytest.fixture(autouse=True, scope="function")
def unique_xdist_group(request):
    # Unique group per test function (nodeid is unique, e.g., test_file.py::test_func)
    group_name = f"boutpp_isolated_{request.node.nodeid.replace('/', '_').replace('::', '_')}"
    request.node.add_marker(pytest.mark.xdist_group(name=group_name))

@pytest.fixture(autouse=True, scope="function")
def cleanup_boutpp():
    yield
    import boutpp
    boutpp.finalise()
