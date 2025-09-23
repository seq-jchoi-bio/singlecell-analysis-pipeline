import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

# -------------------------
# PARAMETERS
# -------------------------
bed_file = sys.argv[1]  # 형식: chr start end name score
output_bed = sys.argv[2]
max_gap = 20
min_cluster_size = 1

# -------------------------
# Load BED and compute summit
# -------------------------
bed = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'name', 'score'])
bed['summit'] = bed['end']  

# -------------------------
# Clustering function
# -------------------------
def cluster_summits(summits, max_gap=20):
    clusters = []
    current = [summits[0]]
    for i in range(1, len(summits)):
        if summits[i] - summits[i-1] <= max_gap:
            current.append(summits[i])
        else:
            clusters.append(current)
            current = [summits[i]]
    clusters.append(current)
    return clusters

# -------------------------
# MM centroid calculation
# -------------------------
def mm_centroid_distance_weighted(summits, max_iter=100, tol=1e-6):
    c = np.mean(summits)
    for _ in range(max_iter):
        distances = np.abs(summits - c) + 1e-6
        weights = 1 / distances
        new_c = np.average(summits, weights=weights)
        if abs(new_c - c) < tol:
            break
        c = new_c
    return c

# -------------------------
# Distance histogram + Cluster count
# -------------------------
all_distances = []
total_cluster_count = 0
valid_cluster_count = 0
bed_records = []

for chrom in bed['chr'].unique():
    chr_bed = bed[bed['chr'] == chrom].sort_values('summit').reset_index(drop=True)
    summits = chr_bed['summit'].values
    clusters = cluster_summits(summits, max_gap=max_gap)

    chrom_total = len(clusters)
    chrom_valid = sum(1 for cl in clusters if len(cl) >= min_cluster_size)
    total_cluster_count += chrom_total
    valid_cluster_count += chrom_valid

    for cl in clusters:
        if len(cl) >= min_cluster_size:
            cl_array = np.array(cl)

            # 중심점 계산
            centroid = mm_centroid_distance_weighted(cl_array)

            # NaN 체크
            if np.isnan(centroid):
                continue

            distances = np.abs(cl_array - centroid)
            all_distances.extend(distances)

            # 점수 계산
            cl_logp = chr_bed[chr_bed['summit'].isin(cl_array)]['score'].values
            score = round(np.sum(cl_logp), 4)

            # BED 레코드 저장
            centroid = int(round(centroid))
            bed_records.append([chrom, centroid, centroid + 1, 'MM_centroid_20', score])

print(f"[INFO] 전체 클러스터 수: {total_cluster_count}")
print(f"[INFO] 조건 충족 (summit ≥ {min_cluster_size}) 클러스터 수: {valid_cluster_count}")

# -------------------------
# Plot histogram
# -------------------------
plt.figure(figsize=(10,6))
plt.hist(all_distances, bins=50, color='steelblue', edgecolor='black')
plt.title("Distribution of distances to MM centroids")
plt.xlabel("Distance to cluster centroid (bp)")
plt.ylabel("Number of summits")
plt.tight_layout()
plt.savefig("MM_centroid_20_distance_histogram.pdf")
plt.show()

# -------------------------
# Save centroid BED file
# -------------------------
bed_df = pd.DataFrame(bed_records, columns=['chr', 'start', 'end', 'name', 'score'])
bed_df.to_csv(output_bed, sep='\t', header=False, index=False)