"""MCP frame-validation tests: a line that parses as JSON but is not a
valid JSON-RPC request must still get an error response for its id — a
client blocked reading that id would otherwise hang forever. Silent-ignore
remains correct only for non-JSON lines and id-less messages.
"""

import json
from io import StringIO
from typing import Any

from gapit.mcp import serve


def exchange(*lines_or_messages: object) -> list[dict[str, Any]]:
    """Feed the loop one line per item (str = raw line, dict = JSON-RPC
    message); return the parsed response objects (test_mcp.py recipe)."""
    raw = "".join(
        (item if isinstance(item, str) else json.dumps(item)) + "\n" for item in lines_or_messages
    )
    out = StringIO()
    serve(StringIO(raw), out)
    return [json.loads(line) for line in out.getvalue().splitlines()]


def test_array_params_gets_32602_not_silence() -> None:
    """Given tools/call whose params is a JSON array (parses, fails frame
    validation), When served, Then a -32602 error answers that id."""
    (response,) = exchange(
        {"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": ["screen"]}
    )
    assert response["id"] == 3
    assert response["error"]["code"] == -32602


def test_request_without_method_gets_32600_not_silence() -> None:
    """Given an object with an id but no method (a valid frame with nothing
    to dispatch), When served, Then a -32600 error answers that id."""
    (response,) = exchange({"jsonrpc": "2.0", "id": 1})
    assert response["id"] == 1
    assert response["error"]["code"] == -32600


def test_non_string_method_gets_32600() -> None:
    """Given a request whose method is a number, When served, Then a -32600
    error answers that id."""
    (response,) = exchange({"jsonrpc": "2.0", "id": 4, "method": 42})
    assert response["id"] == 4
    assert response["error"]["code"] == -32600


def test_invalid_idless_messages_stay_ignored() -> None:
    """Given invalid-shape messages WITHOUT an id around one valid request,
    When served, Then they are ignored like notifications and only the valid
    request is answered."""
    responses = exchange(
        {"method": 42},
        {"jsonrpc": "2.0", "method": "tools/call", "params": ["screen"]},
        {"jsonrpc": "2.0", "id": 9, "method": "tools/list"},
    )
    assert [response["id"] for response in responses] == [9]
