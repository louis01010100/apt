from apt import config


def test_apt_status_line():

    fixture = "#    info     1 816 | 0 error(s) and 670791 warning(s)."
    match = config.APT_STATUS_LINE.match(fixture)

    assert match
    assert 0 == int(match.group(1))
