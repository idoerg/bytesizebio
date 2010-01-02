import UserDict
import MySQLdb
# (c) 2010 Iddo Friedberg
# Distributed under the Biopython license
# http://www.biopython.org/DIST/LICENSE
class GOA_ga(UserDict.UserDict):
    def __init__(self):
        self.data = {}

    def godb_open(self):
        self.dbcon = MySQLdb.connect(host="mysql.ebi.ac.uk", user="go_select",
                   passwd="amigo", db="go_latest", port=4085)
        self.dbcursor = self.dbcon.cursor()

    def godb_close(self):
        self.dbcursor.close()

    def read_gene_assoc(self,inpath):
        for inline in file(inpath):
            db, db_object_id, db_object_symbol, qualifier, go_id, \
            db_reference, evidence, withit, aspect, \
            db_object_name, synonym, db_object_type, \
            taxon_id, date, assigned_by = inline.strip().split('\t')
            key = db+":"+db_object_id
            self.data.setdefault(key,[]).append({'db':db,'db_object_id':db_object_id,
                                            'db_object_symbol': db_object_symbol,
                                            'qualifier': qualifier,
                                            'go_id': go_id,
                                            'db_reference': db_reference,
                                            'evidence': evidence,
                                            'with': withit,
                                            'aspect': aspect,
                                            'db_object_name': db_object_name,
                                            'synonym': synonym,
                                            'db_object_type': db_object_type,
                                            'taxon_id': taxon_id,
                                            'date': date,
                                            'assigned_by': assigned_by})


def find_go_descendants(goa_ga, parent_go_id):
   find_descendants_go = """
   SELECT
     rchild.acc
   FROM
     term AS rchild, term AS ancestor, graph_path
   WHERE
     graph_path.term2_id = rchild.id and
     graph_path.term1_id = ancestor.id and
     ancestor.acc = '%s'
   """
   go_children = {}
   found_genes = {}
   # query GO database for all the child terms of the parent_go_id
   goa_ga.dbcursor.execute(find_descendants_go % parent_go_id)
   rpl = goa_ga.dbcursor.fetchall()
   # Flatten the list. Put it in dictionary keys for faster searching.
   for i in rpl:
      go_children[i[0]] = None
   for gene_id in goa_ga.keys():
      for goa_rec in goa_ga[gene_id]:
         if goa_rec['go_id'] in go_children.keys():
            found_genes.setdefault(gene_id,[]).append(goa_rec)
   return found_genes
