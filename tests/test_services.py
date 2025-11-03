import json
from src.services.metadata_service import MetadataService
from src.services.result_builder import build_result_card


def test_metadata_service_preload_and_reload(tmp_path, monkeypatch):
    data = {"1": {"canonical_3d_smiles": "N/A", "canonical_2d_smiles": "C"}}
    p = tmp_path / "metadata.json"
    p.write_text(json.dumps(data))

    svc = MetadataService(p.as_posix())
    svc.preload()
    assert svc.get_smiles(1) == "C"

    # Update file and reload
    p.write_text(json.dumps({"1": {"canonical_3d_smiles": "CC", "canonical_2d_smiles": "C"}}))
    svc.reload()
    assert svc.get_smiles(1) == "CC"


def test_build_result_card_links():
    entry = {
        'canonical_3d_smiles': 'C',
        'coconut': {'coconut_id': 'COCO123', 'name': 'Foo'},
        'lotus': {'lotus_id': 'LOTUS123', 'name': 'Lotus'},
        'npmrd': {'npmrd_id': 'NP123', 'name': 'NP'}
    }
    card = build_result_card(5, entry, 0.9, None)
    assert card['index'] == 5
    assert card['smiles'] == 'C'
    assert 'coconut' in card['database_links']

