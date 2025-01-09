from src.add_combined_node_attributes import combine_charge_element_to_node
from synutility.SynIO.data_type import load_from_pickle
from src.rc_extract import get_rc_updated
import pytest


def test_combined_charge_element():
    expected_result = {35: "H, 0", 11: "N, 0", 28: "C, 0", 29: "Br, 0"}
    data = load_from_pickle("data/ITS_graphs.pkl.gz")[0]["ITS"]
    data_rc_centre = get_rc_updated(data)
    combine_charge_element_to_node(data_rc_centre)
    try:
        # Attempt to access a valid key (element_charge attribute)
        data_rc_centre_charge_element = data_rc_centre.nodes.data("element_charge")
        assert dict(data_rc_centre_charge_element) == expected_result
    except KeyError:
        pytest.fail("KeyError was raised. No charge_element key")
