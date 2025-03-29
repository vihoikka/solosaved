"""
A simplified locus visualiser that uses:
1. A .gff file from Prokka
2. A HMM hit table from hmmscan
"""

from BCBio import GFF
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os
import gffutils

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# create args inputs for cas_operons_file, locus, output_folder, host_genomes_folder, mode, hmm_rows, hmm_targets, catyper_out, cctyper_protein_table
parser = argparse.ArgumentParser(description="Visualises expanded CRISPR loci")
# parser.add_argument('-i', '--cas_operons_file', help='input file', required=True)
# parser.add_argument('-c', '--cctyper_folder', help='cctyper folder', required=True)
# parser.add_argument('-l', '--locus', help='locus', required=True)
parser.add_argument("-o", "--output_folder", help="output folder", required=True)
# parser.add_argument('-hg', '--host_genomes_folder', help='host genomes folder', required=True)
# parser.add_argument('-v', '--validated_effectors', help='validated effectors', required=True)
# parser.add_argument('-cc', '--cctyper_protein_table', help='cctyper plottable table', required=True)
# parser.add_argument('-k', '--known_effector_table', help='known effectors', required=True)
# parser.add_argument('-rn', '--ring_nucleases', help='ring nucleases', required=True)

parser.add_argument("-pg", "--prokka_gff", help="prokka gff", required=True)
parser.add_argument("-hm", "--hmm_hits", help="hmm hits", required=True)
parser.add_argument("-s", "--sample", help="hmm hits", required=True)
# genbank files
parser.add_argument("-g", "--gbk", help="genbank file", required=True)
parser.add_argument(
    "-go", "--gbk_out_folder", help="genbank output folder", required=True
)

args = parser.parse_args()

gff = args.prokka_gff
hmm = args.hmm_hits
sample = args.sample
output_folder = args.output_folder
gbk = args.gbk
gbk_out_folder = args.gbk_out_folder

# when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the RN
plot_range = 10000

# load the hmm hits file as pandas df. There is no index column
hmm_hits = pd.read_csv(hmm, header=0, index_col=False, delim_whitespace=True)

print(hmm_hits)


def create_gff_iterator(gff_file):
    """
    Creates a GFF iterator for a given GFF file and contig.
    Note that the if the in_handle is closed, then the script will crash later when the iterator is used. This is why we leave it open
    """
    in_handle = open(gff_file)
    gff_iterator = GFF.parse(in_handle)
    return gff_iterator


def add_features_from_gff(gff_iterator, plot_start, plot_end):
    """
    Adds GraphicFeature objects to a list from a given GFF iterator
    """
    temp_features_list = []
    for rec in gff_iterator:
        for feature in rec.features:
            if (
                int(feature.location.start) >= plot_start
                and int(feature.location.end) <= plot_end
            ):
                # print(feature)
                # print("Feature start coord: " + str(feature.location.start))
                # print("Feature end coord: " + str(feature.location.end))

                # if type is CDS, then add to list
                if feature.type == "repeat_region":
                    # print("Repeat region found")
                    strand = 0
                    label = "Repeat region"
                    graphic_feature = GraphicFeature(
                        start=int(feature.location.start) - plot_start,
                        end=int(feature.location.end) - plot_start,
                        strand=strand,
                        color="#000000",
                        label=label,
                    )
                elif feature.type == "CDS":
                    label = feature.qualifiers.get("product", [None])[0]
                    # print("Label: " + str(label))
                    # if label is not NoneType and contains the string "unknown" or "hypthetical", then make label None
                    if label is not None and (
                        "unknown" in label or "hypothetical" in label
                    ):
                        label = None
                    strand = 1 if feature.strand == "+" else -1
                    # print("Adding feature to list")
                    graphic_feature = GraphicFeature(
                        start=int(feature.location.start) - plot_start,
                        end=int(feature.location.end) - plot_start,
                        strand=strand,
                        color="#ffcccc",
                        label=label,
                    )
                temp_features_list.append(graphic_feature)
            # print("----------------")
    return temp_features_list


# number of RN hits
num_rn_hits = len(hmm_hits)
print("Found " + str(num_rn_hits) + " RN hits in sample " + sample)
unique_hits = hmm_hits["query_name"].unique()
print("Unique hits: " + str(unique_hits))

if int(len(unique_hits)) < len(hmm_hits):
    print("Warning: duplicate hits found. Removing duplicates")
    hmm_hits = hmm_hits.drop_duplicates(subset="query_name")
    print("Left with " + str(len(hmm_hits)) + " RN hits")


hmm_hits["RN_coordinate_start"] = None
hmm_hits["RN_coordinate_end"] = None
hmm_hits["strand"] = None
hmm_hits["locus"] = None
hmm_hits["similarity_group_clinker"] = None

# for each hit in the hmm_hits file, search for the entry in the gff file and extract coordinates
for index, row in hmm_hits.iterrows():
    RN_coordinate_start = 0
    RN_coordinate_end = 0
    hit_name = row["query_name"]
    crn_name = row["target_name"]
    # print("Hit name: " + str(hit_name))
    # find the hit in the gff file
    with open(gff) as in_handle:
        for rec in GFF.parse(in_handle):
            for feature in rec.features:
                if feature.id == hit_name:
                    # print("Found hit in GFF file")
                    # print(("Hit name: " + str(feature.id)))
                    RN_coordinate_start = feature.location.start
                    RN_coordinate_end = feature.location.end
                    # print(feature)
                    # print(feature.location.start)
                    # print(feature.location.end)
    # add the coordinates to the hmm_hits df
    hmm_hits.at[index, "RN_coordinate_start"] = RN_coordinate_start
    hmm_hits.at[index, "RN_coordinate_end"] = RN_coordinate_end
    # add strand information
    print("Feature strand: " + str(feature.strand))
    if feature.strand == 1:
        hmm_hits.at[index, "strand"] = 1
    elif feature.strand == -1:
        hmm_hits.at[index, "strand"] = -1
    elif feature.strand == 0:
        hmm_hits.at[index, "strand"] = 0

    # create similarity group for clinker (i.e. name of RN)
    hmm_hits.at[index, "similarity_group_clinker"] = crn_name.split("#")[0].split("_")[
        0
    ]

length_of_genome = 0
# get genome length from gff
with open(gff) as in_handle:
    for rec in GFF.parse(in_handle):
        length_of_genome = len(rec.seq)
        break

print(hmm_hits)

# find out if the hmm hits are within X bases of each other
locus_range = (
    10000  # if RNs are within this range of each other, they are in the same locus
)
# loop through the RNs and assign them to loci depending on their distances
rn_counter = 0
for index, row in hmm_hits.iterrows():
    if rn_counter == 0:  # if starting loop we assign locus 1
        hmm_hits.at[index, "locus"] = 1
    else:  # otherwise check if RN is within range of any previous RN
        for i in range(0, rn_counter):
            if (
                abs(row["RN_coordinate_start"] - hmm_hits.at[i, "RN_coordinate_start"])
                < locus_range
            ):  # if RN is within range of previous RN, assign same locus
                hmm_hits.at[index, "locus"] = hmm_hits.at[i, "locus"]
                break
            else:  # otherwise create new locus
                hmm_hits.at[index, "locus"] = hmm_hits.at[i, "locus"] + 1
    rn_counter += 1

print(hmm_hits)

# create a list of panda dataframes, each containing the hits for a single locus
locus_dfs = []
for i in range(1, hmm_hits["locus"].max() + 1):
    locus_df = hmm_hits[hmm_hits["locus"] == i]
    locus_dfs.append(locus_df)

print("Number of loci: " + str(len(locus_dfs)))

# for each locus, create a visualisation
for i in range(0, len(locus_dfs)):
    locus_df = locus_dfs[i]
    print("Creating viz for locus " + str(i))
    features_list = []  # this features list will contain all RNs within the locus
    print(locus_df)

    # determine locus visualisation start and end points by taking the lowest coordinate of the RNs in the locus and subtracting the plot range
    locus_viz_start = locus_df["RN_coordinate_start"].min() - plot_range
    if locus_viz_start < 0:
        locus_viz_start = 0
    locus_viz_end = locus_df["RN_coordinate_end"].max() + plot_range
    if locus_viz_end > length_of_genome:
        locus_viz_end = length_of_genome

    print("Locus boundaries:")
    print("Start: " + str(locus_viz_start))
    print("End: " + str(locus_viz_end))

    # Add features from GFF
    print("Generating gff based visualisation for locus " + str(i))
    gff_iterator = create_gff_iterator(gff)
    gff_features = add_features_from_gff(
        gff_iterator,
        locus_viz_start,
        locus_viz_end,
    )

    features_list.extend(gff_features)

    for index, row in locus_df.iterrows():  # go through each RN in locus
        gff_iterator = create_gff_iterator(gff)
        RN_start = row["RN_coordinate_start"]
        RN_end = row["RN_coordinate_end"]
        label = row["target_name"].split("#")[0].split("_")[0]
        strand = int(row["strand"])

        print("Adding RN to visualisation")
        RN_graphic_feature = GraphicFeature(
            start=RN_start - locus_viz_start,
            end=RN_end - locus_viz_start,
            strand=strand,
            color="#ebd80c",
            label=label,
        )
        features_list.append(RN_graphic_feature)

    print(features_list)

    # Generate plot using DNA Features Viewer
    record = GraphicRecord(
        sequence_length=locus_viz_end - locus_viz_start, features=features_list
    )
    ax, _ = record.plot(figure_width=10)

    # save
    plt.savefig(output_folder + "/" + sample + "_" + str(index) + ".png")

    # slice the input genbank file to the region of interest using biopython
    # read the genbank file
    record = SeqIO.read(gbk, "genbank")
    # slice the record
    sliced_record = record[locus_viz_start:locus_viz_end]
    # save the sliced record
    SeqIO.write(
        sliced_record,
        gbk_out_folder + "/" + sample + "_" + str(i) + ".gbk",
        "genbank",
    )

    # create a gff file from the sliced gbk record if the ring nuclease is crn4
    if "crn4" in crn_name:
        print("Crn4 hit found. Creating gff file around the locus")
        gff_out = gbk_out_folder + "/" + sample + "_" + str(i) + ".gff"
        gff_out_handle = open(gff_out, "w")
        GFF.write([sliced_record], gff_out_handle)
        gff_out_handle.close()

        # also write the amino acid sequence of the crn4 protein to a fasta file by extracting it from the gbk file
        for feature in sliced_record.features:
            if feature.type == "CDS":
                if "crn4" in feature.qualifiers["product"][0]:
                    print(
                        "Writing Crn4 to fasta file from plasmid "
                        + sample
                        + " locus "
                        + str(i)
                    )
                    protein_seq = feature.qualifiers["translation"][0]
                    protein_seq_record = SeqRecord(
                        Seq(protein_seq),
                        id=sample + "_" + str(i) + "_crn4",
                        description="",
                    )
                    SeqIO.write(
                        protein_seq_record,
                        gbk_out_folder + "/" + sample + "_" + str(i) + "_crn4.fasta",
                        "fasta",
                    )


# create mapping file with each query_name and the corresponding similarity_group_clinker
mapping_df = hmm_hits[["query_name", "similarity_group_clinker"]]
# save as csv without header
mapping_df.to_csv(
    output_folder + "/" + sample + "_similarity_groups.csv", index=False, header=False
)

# features_list = []


# gff_iterator = create_gff_iterator(gff)

# # for every RN, create visualisation of genomic neighbourhood
# for index, row in hmm_hits.iterrows():

#     features_list = []

#     RN_start = row["RN_coordinate_start"]
#     RN_end = row["RN_coordinate_end"]
#     strand = row["strand"]
#     label = row["target_name"].split("#")[0]

#     viz_start = RN_start - plot_range
#     if viz_start < 0:
#         viz_start = 0
#     viz_end = RN_end + plot_range
#     if viz_end > length_of_genome:
#         viz_end = length_of_genome
#     strand = 1 if strand == "+" else -1

#     print("Searching for features in range of RN " + str(index))
#     print("Start coordinate: " + str(RN_start - plot_range))
#     print("End coordinate: " + str(RN_end + plot_range))

#     gff_iterator = create_gff_iterator(gff)

#     # Add features from GFF
#     gff_features = add_features_from_gff(
#         gff_iterator,
#         viz_start,
#         viz_end,
#     )

#     features_list.extend(gff_features)

#     RN_graphic_feature = GraphicFeature(
#         start=RN_start - viz_start,
#         end=RN_end - viz_start,
#         strand=strand,
#         color="#ebd80c",
#         label=label,
#     )
#     features_list.append(RN_graphic_feature)

#     print("Features list:")
#     print(features_list)

#     # Generate plot using DNA Features Viewer
#     record = GraphicRecord(sequence_length=viz_end - viz_start, features=features_list)
#     ax, _ = record.plot(figure_width=10)

#     # save
#     plt.savefig(output_folder + "/" + sample + "_" + str(index) + ".png")

#     # slice the input genbank file to the region of interest using biopython
#     # read the genbank file
#     record = SeqIO.read(gbk, "genbank")
#     # slice the record
#     sliced_record = record[viz_start:viz_end]
#     # save the sliced record
#     SeqIO.write(
#         sliced_record,
#         gbk_out_folder + "/" + sample + "_" + str(index) + ".gbk",
#         "genbank",
#     )
