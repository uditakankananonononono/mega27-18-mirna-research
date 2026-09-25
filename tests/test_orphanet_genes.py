"""Hermetic checks on the committed Orphanet test and its XML parser."""
import json, os, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "orphanet_genes.json")))

X = """<JDBOR><DisorderList><Disorder id="1"><DisorderGeneAssociationList>
<DisorderGeneAssociation><Gene><Symbol>AAA</Symbol></Gene><DisorderGeneAssociationType id="1"><Name lang="en">Disease-causing germline mutation(s) (loss of function) in</Name></DisorderGeneAssociationType><DisorderGeneAssociationStatus id="17991"><Name lang="en">Assessed</Name></DisorderGeneAssociationStatus></DisorderGeneAssociation>
<DisorderGeneAssociation><Gene><Symbol>BBB</Symbol></Gene><DisorderGeneAssociationType id="2"><Name lang="en">Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType><DisorderGeneAssociationStatus id="17991"><Name lang="en">Assessed</Name></DisorderGeneAssociationStatus></DisorderGeneAssociation>
<DisorderGeneAssociation><Gene><Symbol>CCC</Symbol></Gene><DisorderGeneAssociationType id="2"><Name lang="en">Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType><DisorderGeneAssociationStatus id="17997"><Name lang="en">Not yet assessed</Name></DisorderGeneAssociationStatus></DisorderGeneAssociation>
<DisorderGeneAssociation><Gene><Symbol>DDD</Symbol></Gene><DisorderGeneAssociationType id="3"><Name lang="en">Candidate gene tested in</Name></DisorderGeneAssociationType><DisorderGeneAssociationStatus id="17991"><Name lang="en">Assessed</Name></DisorderGeneAssociationStatus></DisorderGeneAssociation>
</DisorderGeneAssociationList></Disorder></DisorderList></JDBOR>"""


def test_parser_filters_status_and_type():
    from orphanet_genes import load_orphanet
    p = os.path.join(tempfile.mkdtemp(), "o.xml"); open(p, "w").write(X)
    dis, lof = load_orphanet(p)
    assert dis == {"AAA", "BBB"} and lof == {"AAA"}


def test_verdicts():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"]
    assert J["O1_dis_vs_expr"]["verdict"] == "CONFIRMED" and J["O2_dis_vs_pub"]["verdict"] == "FALSIFIED"
