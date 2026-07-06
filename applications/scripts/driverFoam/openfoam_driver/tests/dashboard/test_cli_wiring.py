from unittest import mock

from openfoam_driver import cli


def test_dashboard_action_accepted_without_entry():
    # --action dashboard must not require --entry; it should call the launcher.
    with mock.patch("openfoam_driver.dashboard.launch.serve") as serve:
        serve.return_value = 0
        rc = cli.main(["dashboard", "--root", "tutorials",
                       "--port", "9999", "--no-open"])
    assert rc == 0
    assert serve.called
    kwargs = serve.call_args.kwargs
    assert kwargs["port"] == 9999
    assert str(kwargs["root"]).endswith("tutorials")
