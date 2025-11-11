import os
import snapatac2 as snap
import pandas as pd
# 프로젝트 디렉터리 설정
proj_dir = "/home/sohyeong/projects/better_atac_pipeline"

# 필터링된 샘플 파일들이 저장된 디렉터리 설정
filtered_dir = os.path.join(proj_dir, "output_so/filtered_samples/")

# SnapATAC2 데이터가 저장된 경로
snap_data_dir = os.path.join(proj_dir, "output_so/cell_qc/doublets/")

for i in range(1, 10):  # sample1부터 sample9까지 반복
    sample = f"sample{i}"
    
    # 필터링된 파일 경로 설정
    filtered_file = os.path.join(filtered_dir, f"{sample}_filtered.csv")
    try:
        filtered_cells_df = pd.read_csv(filtered_file, index_col="barcode")
    except FileNotFoundError:
        print(f"Filtered file not found for {sample}, skipping...")
        continue
    
    # Fragment 파일 경로 설정
    fragment_file = os.path.join(proj_dir, f"/home/sohyeong/projects/better_atac_pipeline/output_so/cell_ranger_output/{sample}/outs/fragments.tsv.gz")
    try:
        data = snap.pp.import_data(
            fragment_file,
            chrom_sizes=snap.genome.hg38,
            file=os.path.join(snap_data_dir, f"{sample}_doublets.h5ad"),  # 저장 파일명에 샘플 번호 추가
            sorted_by_barcode=False,
        )
    except Exception as e:
        print(f"Error importing data for {sample}: {e}")
        continue

    # 필터링된 바코드 적용
    filtered_barcodes = set(filtered_cells_df.index)
    boolmat_filtered = [barcode in filtered_barcodes for barcode in data.obs_names]
    data.subset(obs_indices=boolmat_filtered)

    # 추가 전처리 단계
    snap.pp.add_tile_matrix(data)
    snap.pp.select_features(data, n_features=250000)
    snap.pp.scrublet(data)
    snap.pp.filter_doublets(data)

    # 결과 출력
    print(f"Number of cells after doublet removal for {sample}: {data.shape[0]}")
