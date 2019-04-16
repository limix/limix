import sys

import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_equal,
    assert_equal,
    assert_string_equal,
)

import limix
from limix.io._fetch import fetch

_samples = [
    "HG00111",
    "HG00112",
    "HG00116",
    "HG00121",
    "HG00133",
    "HG00135",
    "HG00142",
    "HG00143",
    "HG00151",
    "HG00152",
    "HG00154",
    "HG00159",
    "HG00160",
    "HG00171",
    "HG00173",
    "HG00179",
    "HG00189",
    "HG00190",
    "HG00232",
    "HG00233",
    "HG00239",
    "HG00245",
    "HG00253",
    "HG00257",
    "HG00263",
    "HG00274",
    "HG00281",
    "HG00284",
    "HG00309",
    "HG00318",
    "HG00319",
    "HG00330",
    "HG00331",
    "HG00332",
    "HG00343",
    "HG00344",
    "HG00351",
    "HG00357",
    "HG00369",
    "HG00382",
    "HG00419",
    "HG00445",
    "HG00463",
    "HG00472",
    "HG00475",
    "HG00476",
    "HG00531",
    "HG00536",
    "HG00542",
    "HG00553",
    "HG00566",
    "HG00593",
    "HG00595",
    "HG00607",
    "HG00614",
    "HG00625",
    "HG00634",
    "HG00641",
    "HG00650",
    "HG00662",
    "HG00690",
    "HG00702",
    "HG00740",
    "HG01051",
    "HG01055",
    "HG01060",
    "HG01061",
    "HG01069",
    "HG01070",
    "HG01082",
    "HG01098",
    "HG01102",
    "HG01105",
    "HG01113",
    "HG01148",
    "HG01149",
    "HG01171",
    "HG01176",
    "HG01187",
    "HG01188",
    "HG01190",
    "HG01191",
    "HG01204",
    "HG01251",
    "HG01257",
    "HG01259",
    "HG01272",
    "HG01274",
    "HG01344",
    "HG01345",
    "HG01351",
    "HG01357",
    "HG01374",
    "HG01377",
    "HG01378",
    "HG01437",
    "HG01441",
    "HG01492",
    "HG01498",
    "HG01515",
    "HG01516",
    "HG01518",
    "HG01519",
    "HG01521",
    "HG01522",
    "HG01551",
    "HG01617",
    "HG01618",
    "HG01619",
    "HG01620",
    "HG01623",
    "HG01624",
    "HG01625",
    "HG01626",
    "NA06984",
    "NA07051",
    "NA10851",
    "NA11894",
    "NA11932",
    "NA11993",
    "NA11994",
    "NA12273",
    "NA12282",
    "NA12341",
    "NA12489",
    "NA12717",
    "NA12748",
    "NA12750",
    "NA12751",
    "NA12763",
    "NA12814",
    "NA12815",
    "NA12827",
    "NA12830",
    "NA18486",
    "NA18498",
    "NA18502",
    "NA18538",
    "NA18559",
    "NA18560",
    "NA18566",
    "NA18567",
    "NA18571",
    "NA18603",
    "NA18609",
    "NA18613",
    "NA18617",
    "NA18626",
    "NA18631",
    "NA18632",
    "NA18633",
    "NA18636",
    "NA18641",
    "NA18643",
    "NA18647",
    "NA18745",
    "NA18748",
    "NA18858",
    "NA18861",
    "NA18867",
    "NA18871",
    "NA18924",
    "NA18934",
    "NA18943",
    "NA18944",
    "NA18960",
    "NA18965",
    "NA18968",
    "NA18982",
    "NA18987",
    "NA18990",
    "NA18992",
    "NA19002",
    "NA19004",
    "NA19020",
    "NA19041",
    "NA19057",
    "NA19060",
    "NA19067",
    "NA19068",
    "NA19075",
    "NA19077",
    "NA19080",
    "NA19082",
    "NA19088",
    "NA19102",
    "NA19116",
    "NA19121",
    "NA19130",
    "NA19146",
    "NA19147",
    "NA19160",
    "NA19172",
    "NA19190",
    "NA19197",
    "NA19204",
    "NA19308",
    "NA19309",
    "NA19310",
    "NA19312",
    "NA19347",
    "NA19351",
    "NA19352",
    "NA19355",
    "NA19360",
    "NA19397",
    "NA19403",
    "NA19428",
    "NA19430",
    "NA19451",
    "NA19453",
    "NA19463",
    "NA19467",
    "NA19472",
    "NA19625",
    "NA19648",
    "NA19654",
    "NA19655",
    "NA19657",
    "NA19660",
    "NA19661",
    "NA19664",
    "NA19717",
    "NA19723",
    "NA19725",
    "NA19731",
    "NA19734",
    "NA19740",
    "NA19749",
    "NA19756",
    "NA19761",
    "NA19771",
    "NA19774",
    "NA19789",
    "NA19795",
    "NA19819",
    "NA19834",
    "NA19908",
    "NA19920",
    "NA20127",
    "NA20278",
    "NA20281",
    "NA20287",
    "NA20291",
    "NA20314",
    "NA20334",
    "NA20339",
    "NA20341",
    "NA20344",
    "NA20348",
    "NA20357",
    "NA20359",
    "NA20412",
    "NA20414",
    "NA20505",
    "NA20507",
    "NA20508",
    "NA20517",
    "NA20518",
    "NA20521",
    "NA20525",
    "NA20527",
    "NA20534",
    "NA20537",
    "NA20581",
    "NA20582",
    "NA20753",
    "NA20754",
    "NA20768",
    "NA20771",
    "NA20772",
    "NA20774",
    "NA20775",
    "NA20804",
]


def test_fetch_csv():

    with limix.file_example("expr.csv") as filepath:

        spec = f"{filepath}:csv:row=trait,trait[gene1]"
        y = fetch("trait", spec, verbose=False)

        assert_string_equal(y.name, "trait")
        assert_array_equal(y["sample"], _samples)
        assert_equal(y.shape, (274, 1))
        assert_allclose(y.values[:2, 0], [-3.752_345_147_31, -0.421_128_991_488])
        assert_array_equal(y.coords, ["trait", "sample"])

        spec = f"{filepath}:csv:row=trait,trait[gene11]"
        y = fetch("trait", spec, verbose=False)

        assert_string_equal(y.name, "trait")
        assert_array_equal(y["sample"], _samples)
        assert_equal(y.shape, (274, 1))
        assert_allclose(y.values[:2, 0], [0.798_312_717_19, 0.237_496_587_19])

        spec = f"{filepath}:csv:row=trait"
        y = fetch("trait", spec, verbose=False)

        assert_string_equal(y.name, "trait")
        assert_array_equal(y["sample"], _samples)
        assert_equal(y.shape, (274, 11))

        spec = f"{filepath}:csv:row=trait,col=sample"
        y = fetch("trait", spec, verbose=False)
        assert_equal(y.shape, (274, 11))
        assert_equal(y.dims, ("sample", "trait"))

        spec = f"{filepath}:csv:row=sample,col=trait"
        y = fetch("trait", spec, verbose=False)
        assert_equal(y.shape, (11, 274))
        assert_equal(y.dims, ("sample", "trait"))

        spec = f"{filepath}:csv:"
        y = fetch("trait", spec, verbose=False)

        spec = f"{filepath}:csv"
        y = fetch("trait", spec, verbose=False)

        spec = f"{filepath}:csv:row=samples"
        with pytest.raises(ValueError):
            y = fetch("trait", spec, verbose=False)

        spec = "wrong_filepath:csv:row=sample"
        with pytest.raises(FileNotFoundError):
            y = fetch("trait", spec, verbose=False)

        spec = f"{filepath}:csv:row=sample,col=trait"
        with pytest.raises(ValueError):
            y = fetch("traits", spec, verbose=False)

        spec = f"{filepath}:csvs:row=sample,col=trait"
        with pytest.raises(ValueError):
            y = fetch("trait", spec, verbose=False)


def test_fetch_bed():
    import os

    filenames = [
        "chrom22_subsample20_maf0.10.bed",
        "chrom22_subsample20_maf0.10.fam",
        "chrom22_subsample20_maf0.10.bim",
    ]

    with limix.file_example(filenames) as filepaths:
        folder = os.path.dirname(filepaths[0])
        filepath = os.path.join(folder, "chrom22_subsample20_maf0.10")

        specs = [
            f"{filepath}",
            f"{filepath}:bed:row=candidate,col=sample",
            f"{filepath}:bed:col=sample",
            f"{filepath}:bed",
            f"{filepath}",
            f"{filepath}:",
            f"{filepath}::",
            f"{filepath}:bed:",
            f"{filepath}::row=candidate",
        ]

        for spec in specs:
            G = fetch("genotype", spec, verbose=False)

            assert_string_equal(G.name, "genotype")
            assert_array_equal(G["sample"], _samples)
            assert_equal(G.shape, (274, 49008))
            assert_array_equal(G.dims, ["sample", "candidate"])


@pytest.mark.skipif(sys.platform.startswith("win"), reason="for unix only")
def test_fetch_win_drive_on_unix():
    with pytest.raises(ValueError):
        fetch("trait", r"C:\Temp\Wrong.csv:csv:row=sample", verbose=False)


@pytest.mark.skipif(not sys.platform.startswith("win"), reason="for win only")
def test_fetch_win_drive_on_win():
    with pytest.raises(FileNotFoundError):
        fetch("trait", r"C:\Temp\Wrong.csv:csv:row=sample", verbose=False)
