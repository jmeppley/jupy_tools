from edl import taxon as edltaxon, util
import os, re, time
from Bio import Entrez
from IPython.core.interactiveshell import getoutput
from requests import HTTPError
import logging


class TaxonomicDatabase:
    """ Object representing a sequence db with taxonomic info. EG RefSeq """

    def __init__(self, *db_path):
        """ Loads the db files into RAM 

        Arguments can be a single db prefix of a fully formd py-metagenomics db:

         EG: /mnt/seqdbs/RefSeq/200/RefSeq.ldb/lastdb
         where the files are:
            taxid map: {prefix}.tax
            taxdump dir: dirname(prefix)

        Or two paths:

            taxid_map location
            taxdump dir

        """
        self.db_path = db_path
        if isinstance(db_path, str) or len(db_path) == 1:
            if not isinstance(db_path, str):
                db_path = next(iter(db_path))

            self.db_dir = os.path.dirname(db_path)
            self.db_taxid_map = db_path + ".tax"
        else:
            self.db_taxid_map , self.db_dir = db_path
        self.taxonomy = edltaxon.readTaxonomy(self.db_dir)
        self.taxid_map = util.parseMapFile(self.db_taxid_map, valueType=int)

    def find_missing_taxids(self):
        """
        Inspect the database files to see if there are any taxids in the accession map
        that are missing from the taxonomy.
        """
        taxids_from_accs = set(self.taxid_map.values())
        taxids_from_tax = set(self.taxonomy.idMap.keys())
        return taxids_from_accs.difference(taxids_from_tax)

    def fix_taxid_map(self, new_taxid_map):
        """ modify the taxids in the acc to taxid map """
        self.db_taxid_map_bak = fix_taxid_map(new_taxid_map, self.db_taxid_map)


def _backup_taxid_map(taxid_map_file):
    """ Make a copy of the taxid map """
    n = 0
    taxid_map_file_bak = taxid_map_file + ".{}.bak".format(n)
    while os.path.exists(taxid_map_file_bak):
        n += 1
        taxid_map_file_bak = taxid_map_file + ".{}.bak".format(n)
    message = getoutput(f"cp {taxid_map_file} {taxid_map_file_bak}")
    if len(message) > 0:
        logging.warning(message)
    return taxid_map_file_bak


def fix_taxid_map(new_taxid_map, taxid_map_file):
    """ modify the taxids in the acc to taxid map """
    taxid_map_file_bak = _backup_taxid_map(taxid_map_file)
    if not os.path.exists(taxid_map_file_bak):
        raise Exception("Backup did not work! Aborting.")
    taxid_translate_expressions = get_subst_expressions(new_taxid_map)
    translate_file(taxid_map_file_bak, taxid_map_file, taxid_translate_expressions)
    return taxid_map_file_bak


def get_subst_expressions(word_map, key_on_old=True):
    translate_expressions = []
    for old, new in word_map.items():
        if not key_on_old:
            new, old = (old, new)
        translate_expressions.append(r"s/\b{}\b/{}/".format(old, new))
    return translate_expressions


def translate_file(source, destination, expressions, multiple_pipes=False):
    if multiple_pipes:
        command = " | ".join(["perl -pe'" + ex + "'" for ex in expressions])
    else:
        command = "perl -pe '" + "; ".join(expressions) + "'"

    return getoutput(f"cat {source} | {command} > {destination}")


def lookup_new_taxids(old_taxids, database=None, sleep=0, email="jmeppley@hawaii.edu"):
    """ For a list of deprecated taxids, get the new taxid from Entrez """
    Entrez.email = email

    new_taxid = {}
    count = 0
    query_count = 0
    for taxid in old_taxids:
        count += 1
        try:
            if query_count % 3 == 2:
                time.sleep(sleep)
            query_count += 1
            # first see if the taxid redirects somewhere
            handle = Entrez.esummary(db="taxonomy", id=taxid)
            r = handle.read()
            m = re.search(r"AkaTax.+\>(\d+)\<\/Item", r)
            if m:
                aka_taxid = int(m.group(1))
                if aka_taxid != 0:
                    new_taxid[taxid] = aka_taxid
                    continue

            if database is None:
                continue

            # try looking up an accession, if user gave us a db object
            for acc, acc_taxid in database.taxid_map.items():
                if acc_taxid == taxid:
                    if query_count % 3 == 2:
                        time.sleep(sleep)
                    query_count += 1
                    handle = Entrez.esummary(db="protein", id=acc)
                    r = handle.read()
                    m = re.search(r'Name="TaxId"[^>]+>(\d+)<', r)
                    if m:
                        new_taxid[taxid] = int(m.group(1))
                        break

        except Exception as error:
            print(f"Exception on {count}th taxid")
            raise
    return new_taxid
