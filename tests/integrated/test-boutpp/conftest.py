import pytest

@pytest.fixture(autouse=True, scope="function")
def xdist_each(request):
    request.node.add_marker(pytest.mark.xdist_group(name="each_test"))  # Forces one per worker if combined

