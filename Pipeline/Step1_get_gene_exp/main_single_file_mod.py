import sys
import os
from intermine.webservice import Service




def main():



    if not os.path.exists("results"):
        os.makedirs("results")


    service = Service("https://apps.araport.org/thalemine/service")
    file = open("results/all_genes.csv", "w")
    list_written = []
    list_genes = []
    list_gene_names = []
    for index , line in enumerate(open(os.getcwd()+"/" + sys.argv[1])):
        gene = line.strip()
        query = service.new_query("Gene")
        query.add_view(
            "primaryIdentifier", "RNASeqExpressions.expressionLevel",
            "RNASeqExpressions.experiment.SRAaccession",
            "RNASeqExpressions.experiment.tissue", "RNASeqExpressions.unit"
        )
        query.add_sort_order("Gene.RNASeqExpressions.experiment.SRAaccession", "DESC")
        query.add_constraint("primaryIdentifier", "=", gene, code="A")



        for row in query.rows():
            experiment_tissue = str(row["RNASeqExpressions.experiment.SRAaccession"]) + "-" + str(row["RNASeqExpressions.experiment.tissue"])
            expression_value = str(row["RNASeqExpressions.expressionLevel"])
            if experiment_tissue not in list_written:
                list_written.append(experiment_tissue)
            list_genes.append((gene , experiment_tissue , expression_value))


        list_gene_names.append(gene)




    for item in list_written:
        file.write("\t" + item)

    file.write("\n")


    flag = 0
    for gene_name in list_gene_names:
        file.write(gene_name)
        for item in list_written:
            flag = 0
            for gene_name_temp , exp_tissue , expression_value in list_genes:
                if gene_name == gene_name_temp:
                    if item == exp_tissue:
                        file.write("\t" + expression_value)
                        flag = 1
                        break
            if flag == 0:
                file.write("\t0")

        file.write("\n")



    file.close()

















if __name__ == '__main__':
   main()