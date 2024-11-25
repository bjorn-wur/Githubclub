from bif_template.main import return_hello_world


def test_return_hello_world():
    assert return_hello_world() == "Hello World!"
