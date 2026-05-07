from __future__ import print_function
import os
import sys

input_file = r"C:/毕业论文/数据处理及绘图/ppi_build/9606.protein.links.full.v12.txt"
output_file = r"C:/毕业论文/数据处理及绘图/ppi_build/PPIN.txt"

if not os.path.exists(input_file):
    sys.exit("Can't locate input file %s" % input_file)

prior = 0.041

def compute_prior_away(score, prior):
    if score < prior:
        score = prior
    score_no_prior = (score - prior) / (1 - prior)
    return score_no_prior

header = True

#  打开输出文件
with open(output_file, "w") as out:

    # 可选：写表头（强烈推荐）
    out.write("protein1\tprotein2\tcombined_score\n")

    for line in open(input_file):

        if header:
            header = False
            continue
        
        l = line.split()
        
        (protein1, protein2,
         neighborhood, neighborhood_transferred,
         fusion, cooccurrence,
         homology,
         coexpression, coexpression_transferred,
         experiments, experiments_transferred,
         database, database_transferred,
         textmining, textmining_transferred,
         initial_combined) = l

        ## divide by 1000
        neighborhood = float(neighborhood) / 1000
        neighborhood_transferred = float(neighborhood_transferred) / 1000
        fusion = float(fusion) / 1000
        cooccurrence = float(cooccurrence) / 1000
        homology = float(homology) / 1000
        coexpression = float(coexpression) / 1000
        coexpression_transferred = float(coexpression_transferred) / 1000
        experiments = float(experiments) / 1000
        experiments_transferred = float(experiments_transferred) / 1000
        database = float(database) / 1000
        database_transferred = float(database_transferred) / 1000
        #textmining_transferred = float(textmining_transferred) / 1000
        initial_combined = int(initial_combined)

        ## compute prior away
        neighborhood_prior_corrected = compute_prior_away(neighborhood, prior)
        neighborhood_transferred_prior_corrected = compute_prior_away(neighborhood_transferred, prior)
        fusion_prior_corrected = compute_prior_away(fusion, prior)
        cooccurrence_prior_corrected = compute_prior_away(cooccurrence, prior)
        coexpression_prior_corrected = compute_prior_away(coexpression, prior)
        coexpression_transferred_prior_corrected = compute_prior_away(coexpression_transferred, prior)
        experiments_prior_corrected = compute_prior_away(experiments, prior)
        experiments_transferred_prior_corrected = compute_prior_away(experiments_transferred, prior)
        database_prior_corrected = compute_prior_away(database, prior)
        database_transferred_prior_corrected = compute_prior_away(database_transferred, prior)
        #textmining_prior_corrected                   = compute_prior_away (textmining, prior)            
        #textmining_transferred_prior_corrected       = compute_prior_away (textmining_transferred, prior) 

        ## combine direct + transferred
        neighborhood_both_prior_corrected = 1.0 - (1.0 - neighborhood_prior_corrected) * (1.0 - neighborhood_transferred_prior_corrected)
        coexpression_both_prior_corrected = 1.0 - (1.0 - coexpression_prior_corrected) * (1.0 - coexpression_transferred_prior_corrected)
        experiments_both_prior_corrected = 1.0 - (1.0 - experiments_prior_corrected) * (1.0 - experiments_transferred_prior_corrected)
        database_both_prior_corrected = 1.0 - (1.0 - database_prior_corrected) * (1.0 - database_transferred_prior_corrected)

        ## combine scores (exclude textmining)
        combined_score_one_minus = (
            (1.0 - neighborhood_both_prior_corrected) *
            (1.0 - fusion_prior_corrected) *
            (1.0 - cooccurrence_prior_corrected) *
            (1.0 - coexpression_both_prior_corrected) *
            (1.0 - experiments_both_prior_corrected) *
            (1.0 - database_both_prior_corrected)
        )

        combined_score = (1.0 - combined_score_one_minus)
        combined_score *= (1.0 - prior)
        combined_score += prior

        combined_score = int(combined_score * 1000)

        #  写入文件（制表符分隔）
        out.write(f"{protein1}\t{protein2}\t{combined_score}\n")

print("Finished. Output written to:", output_file)
