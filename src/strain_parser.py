#!/usr/bin/env python3

"""Parses FluClust results for given strains.

Usage:
  strain_parser.py --fasta -q <strains> -c <clusters> -o <output>
  strain_parser.py --text -q <strains> -c <clusters> -o <output>

Options:
  -h --help                         Show this screen.
  --fasta                           Run the script in fasta mode.
  --text                            Run the script in text mode.
  -q --query=<strains>              Comma-separated list of strains to parse.
                                    In fasta mode the query is a fasta file that
                                    is blasted agains the FluClust results.
                                    In text mode the query is a text file with
                                    names of the strains to be queried agains
                                    the FluClust results.
  -c --clusters=<clusters>          FluClust clusters filepath.
  -o --output=<output>              Output filepath.
"""


from docopt import docopt
import glob
import thefuzz


def query_text(query, clusters, fuzy_search=False):
    """Query the clusters file for the given strains.

    Args:
        query (str): Comma-separated list of strains to query.
        clusters (str): Path to the clusters file.

    Returns:
        list: List of the strains that were found in the clusters file.
    """
    hits = {}
    cluster_ids = glob.glob(f"{clusters}/*/cluster_influenza_segment_*.csv")
    cluster_metadata = glob.glob(f"{clusters}/*/meta_influenza_segment_*.csv")

    with open(query) as query_file:
        for query_strain in query_file:
            for cluster in cluster_metadata:
                with open(cluster) as cluster_file:
                    for strain in cluster_file:
                        if strain.strip() in query_strain.strip():
                            hits[query_strain.strip()] = [strain.strip(), cluster]
    return hits
    print(hits)


def main():
    arguments = docopt(__doc__)
    if arguments["--fasta"]:
        pass
    if arguments["--text"]:
        query_text(arguments["--query"], arguments["--clusters"])


if __name__ == "__main__":
    main()
