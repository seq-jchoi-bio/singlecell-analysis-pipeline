# ------------------------------------------------------------------------------------
# [1] 더블릿 제거된 H5AD 파일 이름 변경 (안전하게 '복사'로 수행)
#  - 원본 파일을 건드리지 않도록 shutil.copy를 사용하여 복사본 생성
#  - sample1 ~ sample9 까지 순회하며,
#    /doublets/{sample}_doublets.h5ad → /aggr/{sample}.h5ad 로 복사.
# ------------------------------------------------------------------------------------
#Change Doublet removal H5AD files' names
#파일 속 내용 깨질 수 있으므로 안전하게 복사본 사용
import shutil

for i in range(1, 10):  
    sample = f"sample{i}"  
    # 복사 작업 수행 (원본은 그대로 두고, 대상 디렉터리에 동일 내용으로 새 파일 생성)
    shutil.copy(
        f"/home/sohyeong/projects/better_atac_pipeline/output_so/cell_qc/doublets/{sample}_doublets.h5ad",  
        f"/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/{sample}.h5ad" 
    )
    print(f"{sample} 파일 복사 완료")  


# ------------------------------------------------------------------------------------
# [2] H5AD 파일들을 하나의 AnnDataSet(.h5ads)으로 병합
#  - snapatac2의 AnnDataSet은 다수의 per-sample h5ad를 묶어 관리하기 좋음.
#  - adatas 인자: (그룹명, 파일경로) 튜플 목록
#  - 병합 결과는 단일 .h5ads 파일에 저장.
# ------------------------------------------------------------------------------------
# Merge all H5AD files
# 개별 H5AD 파일들을 하나의 파일로 병합
import snapatac2 as snap

# 병합에 사용할 파일 경로 목록 (sample1.h5ad ~ sample9.h5ad)
file_paths = [
    f"/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/sample{i}.h5ad"
    for i in range(1, 10)
]

# 병합된 데이터 저장 경로
output_path = "/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/colon.h5ads"

# 병합 및 저장 시도
print("Merging and saving data...")
try:
    data = snap.AnnDataSet(
        adatas=[(f"sample{i+1}", path) for i, path in enumerate(file_paths)],
        filename=output_path
    )
    print(f"Merged data saved to {output_path}")
except Exception as e:
    print(f"Error during data merging: {e}")


# ------------------------------------------------------------------------------------
# [3] 병합된 .h5ads 파일이 정상적으로 로드되는지 확인
#  - snap.read 로 읽고, 기본 메타 정보(셀 수, 피처 수)를 확인
# ------------------------------------------------------------------------------------
#병합된 파일 속 내용 확인
import snapatac2 as snap

# 병합된 데이터 파일 경로
merged_file_path = "/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/colon.h5ads"

# 데이터 로드 및 확인
try:
    print(f"Loading merged dataset from {merged_file_path}...")
    merged_data = snap.read(merged_file_path)
    print("Merged dataset successfully loaded.")
    print(f"Number of cells: {merged_data.n_obs}")
    print(f"Number of features: {merged_data.n_vars}")
except Exception as e:
    print(f"Error loading merged dataset: {e}")


# ------------------------------------------------------------------------------------
# [4] 클러스터링 전처리 및 임베딩
#  - data: 위에서 생성한 AnnDataSet (여기서 data 객체를 계속 사용)
#  - 고변이 피처 선택 → 스펙트럴 임베딩 → 배치 교정(Harmony) → KNN/Leiden → UMAP
#  - 스레드(병렬화) 환경변수: 과도한 스레드 생성으로 OOM/과열 방지 목적
# ------------------------------------------------------------------------------------
#Clustering
data

# 고변이 피처 선택 (ATAC 특성상 피처 수가 많으므로 n_features를 크게)
snap.pp.select_features(data, n_features=200000)

# 스펙트럴 임베딩 계산 (그래프 기반 저차원 표현 학습)
snap.tl.spectral(data)

# 병렬 스레드 제한 설정 (과도한 스레드 생성으로 인한 성능 저하/메모리 폭주 방지)
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"]      = "1"
os.environ["MKL_NUM_THREADS"]      = "1"
os.environ["NUMEXPR_NUM_THREADS"]  = "1"
os.environ["RAYON_NUM_THREADS"]    = "1"  # (snapATAC2의 Rust 병렬화)
os.environ["MALLOC_ARENA_MAX"]     = "2"  # glibc arena 폭주 방지

# 배치 효과 보정 (샘플 간 시스템적 차이를 제거)
#  - batch="sample" 메타데이터 컬럼을 기준으로 Harmony 수행
#  - use_dims=30: 스펙트럴 임베딩의 처음 30개 차원을 사용
snap.pp.harmony(data, batch="sample", use_dims=30, max_iter_harmony=20)

# KNN 그래프 구성 (Harmony로 교정된 임베딩 사용)
snap.pp.knn(data, n_neighbors=10, use_dims=30, use_rep="X_spectral_harmony") # Using embeding from harmony

# Leiden 클러스터링 (재현 가능성을 위해 random_state 고정)
snap.tl.leiden(data, random_state=0)  # Leiden clustering
print("Leiden clustering completed.")

# UMAP 계산 (시각화를 위한 2D 임베딩)
snap.tl.umap(data, use_rep="X_spectral_harmony", use_dims=30)

# UMAP 시각화 (클러스터별 색)
sc_plt = snap.pl.umap(data, interactive=False, color="leiden")
import matplotlib.pyplot as plt
plt.savefig("/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/umap_leiden.png", dpi=300)
plt.savefig("/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/umap_leiden.svg")

# UMAP 시각화 (샘플별 색)
sc_plt = snap.pl.umap(data, interactive=False, color="sample")
plt.savefig("/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/umap_sample.png", dpi=300)
plt.savefig("/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/umap_sample.svg")


# ------------------------------------------------------------------------------------
# [5] 유전자 매트릭스 생성 및 UMAP 
#  - ATAC 데이터에서 gene activity matrix(유전자 수준) 생성
#  - 시각화 일관성을 위해 기존 UMAP 좌표를 gene_matrix에 복사
# ------------------------------------------------------------------------------------
#Cell Annotation
%%time
gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38)
gene_matrix

# 기존 UMAP 좌표를 gene_matrix로 복사하여 동일 좌표계에서 시각화(같은 UMAP 좌표를 복사하면 세포 클러스터링 결과 (bin/peak 기반)와 세포 타입 annotation (gene 기반)을 같은 2D 평면에서 직접 비교 가능
gene_matrix.obsm["X_umap"] = data.obsm["X_umap"]

sc.pl.umap(gene_matrix, use_raw=False, color="leiden", legend_loc="on data")
sc.pl.umap(gene_matrix, use_raw=False, color="sample")

# ------------------------------------------------------------------------------------
# [6] 마커 유전자 정의
#  - 각 세포 타입을 구분하기 위한 대표 마커 목록
#  - 실제 실험/참조 논문에 맞춰 조정 필요
# ------------------------------------------------------------------------------------
import scanpy as sc
import matplotlib.pyplot as plt

# 마커 유전자 정의
marker_genes = {
    "Atrial Cardiomyocytes": ["NPPA", "MYH6", "MYL7"],
    "Ventricular Cardiomyocytes": ["MYH7", "HEY2", "MYL2"],
    "Fibroblasts": ["DCN"],
    "Endothelial": ["EGFL7", "VWF"],
    "Smooth Muscle": ["GJA4", "TAGLN"],
    "Macrophages": ["CD163", "MS4A6A"],
    "Lymphocytes": ["IL7R", "THEMIS"],
    "Adipocytes": ["ADIPOQ", "CIDEA"],
    "Nervous Cells": ["NRXN3", "GPM6B"],
    "Endocardial-like Cells": ["NRG3", "NPR3"],
    "Myofibroblasts": ["MYH10"],
    "Arterial Smooth Muscle Cells": ["ACTA2", "TAGLN"],
    "Epiblast":["NANOG", "POU5F1", "SOX2", "KLF17", "TDGF1"],
    "Primitive Endoderm":["GATA6", "GATA4", "SOX17", "PDGFRA"],
    "Trophectoderm": ["GATA3", "GATA2", "KRT18", "TEAD3"],
    "Early Inner Cell Mass": ["RSPO3", "ARGFX", "PRDM14", "SOX2"]
}

# ------------------------------------------------------------------------------------
# [7] Dot plot으로 마커 발현 시각화
#  - groupby="leiden": 클러스터별 마커 발현 패턴 확인
#  - standard_scale="var": 유전자(열)별 표준화 → 비교 용이
#  - dot_min/max, cmap, figsize 등 시각화 가독성 설정
# ------------------------------------------------------------------------------------
# Dot Plot 그리기
dot = sc.pl.dotplot(
    gene_matrix,  # 데이터 객체
    var_names=marker_genes,  # 마커 유전자
    groupby="leiden",  # 클러스터별 그룹
    standard_scale="var",  # 변수를 표준화
    dot_min=0.1,  # 최소 점 크기
    dot_max=1,  # 최대 점 크기
    cmap="Blues",  # 컬러맵 설정
    figsize=(12, 10), 
    show=False  # show=False로 설정해야 저장 가능
)

# 저장
plt.savefig("/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/dotplot.png", dpi=300, bbox_inches="tight")
plt.savefig("/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/dotplot.svg", bbox_inches="tight")
plt.close()


# ------------------------------------------------------------------------------------
# [8] per-sample h5ad에 클러스터/세포타입 메타데이터 추가 후 업데이트 저장
#  - gene_matrix에서 클러스터(Leiden), 샘플, 바코드 추출
#  - 클러스터→세포타입 매핑(dict) 적용
#  - 각 sampleX.h5ad를 읽어 obs에 'celltype', 'leiden', 'sample', 'sample_cluster' 기록
#  - *_updated.h5ad 로 저장 후, 다시 AnnDataSet(.h5ads)로 묶어 저장
# ------------------------------------------------------------------------------------
#Update H5AD files(Add cluster information)
import os
import snapatac2 as snap
import scanpy as sc
import pandas as pd

# 원본 per-sample h5ad들 (수정 대상)
file_paths = [
    f"/home/sohyeong/projects/better_atac_pipeline/output_so/aggr/sample{i}.h5ad"
    for i in range(1, 10)
]

# 업데이트 저장 경로 (없으면 생성)
aggr_update_dir = "/home/sohyeong/projects/better_atac_pipeline/output_so/aggr_update"
os.makedirs(aggr_update_dir, exist_ok=True)

# gene_matrix: 바코드 중복 제거 (중복 인덱스로 인한 reindex 오류 방지)
gene_matrix = gene_matrix[~gene_matrix.obs.index.duplicated(keep='first')].copy()

# gene_matrix 메타데이터 추출 (문자형으로 통일하여 카테고리 변환 시 혼란 방지)
leiden = gene_matrix.obs["leiden"].astype(str)
sample = gene_matrix.obs["sample"].astype(str)
cell_barcodes = gene_matrix.obs.index.astype(str)

# 클러스터 → 세포 유형 매핑 (알 수 없으면 "Unknown")
cluster_ctype_dic = {
    '2': "Atrial_Cardiomyocytes", '16': "Atrial_Cardiomyocytes", '20': "Atrial_Cardiomyocytes",
    '0': "Ventricular_Cardiomyocytes", '4': "Ventricular_Cardiomyocytes",
    '3': "Ventricular_Cardiomyocytes", '5': "Ventricular_Cardiomyocytes",
    '6': "Ventricular_Cardiomyocytes", '7': "Ventricular_Cardiomyocytes",
    '8': "Ventricular_Cardiomyocytes", '11': "Ventricular_Cardiomyocytes", '15': "Ventricular_Cardiomyocytes",
    '10': "Primitive_Endoderm", '13': "Primitive_Endoderm",
    '12': "Endothelial", '14': "Endothelial",
    '19': "Myofibroblasts", '21': "Nervous_Cells", '22': "Nervous_Cells", '18': "Macrophages", '17': "Smooth_Muscle",
    '1': "Fibroblasts", '9': "Fibroblasts"
}
cell_types = leiden.map(lambda l: cluster_ctype_dic.get(l, "Unknown"))

# gene_matrix 기반 메타 테이블 (인덱스=셀 바코드)
df = pd.DataFrame(
    {"leiden": leiden, "sample": sample, "celltype": cell_types},
    index=cell_barcodes,
)

# per-sample h5ad 업데이트 루프
updated_file_paths = []  # [(group_name, path), ...] → 이후 AnnDataSet 작성에 사용
for h5ad_path in file_paths:
    try:
        print(f"Processing {h5ad_path} ...")
        ad = sc.read_h5ad(h5ad_path)

        # 바코드 중복 제거 (인덱스 유일성 확보)
        ad = ad[~ad.obs.index.duplicated(keep='first')].copy()

        # gene_matrix 메타를 현재 파일 바코드 순서로 정렬 (인덱스 기준 align)
        merged = df.reindex(ad.obs.index.astype(str))

        # 파일명에서 sample 이름 추출 (예: sample3)
        sample_name = os.path.splitext(os.path.basename(h5ad_path))[0]

        # 누락값 보정: 해당 sample에 없는 셀/정보는 기본값으로 채움      
        merged["celltype"] = merged["celltype"].fillna("Unknown").astype(str)
        merged["leiden"]   = merged["leiden"].fillna("NA").astype(str)
        merged["sample"]   = merged["sample"].fillna(sample_name).astype(str)

        # obs 메타데이터 갱신 
        ad.obs["celltype"] = pd.Categorical(merged["celltype"].values)
        ad.obs["leiden"]   = pd.Categorical(merged["leiden"].values)
        ad.obs["sample"]   = pd.Categorical(merged["sample"].values)
        # sample_cluster: 세포타입과 샘플명을 결합한 그룹 레이블 (ex: Ventricular_Cardiomyocytes.sample1)
        ad.obs["sample_cluster"] = ad.obs["celltype"].astype(str) + "." + ad.obs["sample"].astype(str)

        # 저장 (aggr_update/sampleX_updated.h5ad)
        updated_path = os.path.join(aggr_update_dir, f"{sample_name}_updated.h5ad")
        ad.write(updated_path)
        updated_file_paths.append((sample_name, updated_path))
        print(f"✅ Updated: {updated_path} (matched {merged['leiden'].notna().sum()} cells)")

    except Exception as e:
       
        print(f"❌ Error processing {h5ad_path}: {e}")

# 새로운 AnnDataSet(.h5ads) 생성 (업데이트된 per-sample h5ad를 묶기)
updated_dataset_path = os.path.join(aggr_update_dir, "colon_updated.h5ads")
snap.AnnDataSet(adatas=updated_file_paths, filename=updated_dataset_path)
print(f"✅ New AnnDataSet with celltype added: {updated_dataset_path}")

# 사용 완료된 핸들 정리
data.close()



# ------------------------------------------------------------------------------------
# [10] Fragment 파일 내보내기
#  - AnnDataSet(data)와 groupby 시리즈를 기반으로 fragment bed를 샘플×클러스터 단위로 분할 저장
#  - sample_cluster 레이블을 그대로 사용하여 BED 파일 접두/접미사로 기록
# ------------------------------------------------------------------------------------
#Convert to Fragment file
data = snap.AnnDataSet(adatas=updated_file_paths, filename=updated_dataset_path)
print(f"✅ New AnnDataSet with celltype added: {updated_dataset_path}")

import pandas as pd
import scanpy as sc
import snapatac2 as snap
import os

# per-sample에서 sample_cluster 시리즈를 모아 합치기
labels = []
for name, path in updated_file_paths:
    ad = sc.read_h5ad(path)
    if "sample_cluster" not in ad.obs:
        raise ValueError(f"{path}에 sample_cluster가 없습니다.")
    labels.append(ad.obs["sample_cluster"].astype(str))

# 인덱스: 전체 셀 바코드 (각 샘플의 바코드를 그대로 이어붙임)
groupby_series = pd.concat(labels)   # index = 전체 셀 바코드

# 내보내기 디렉터리 준비
dir_fragments = "/home/sohyeong/projects/better_atac_pipeline/output_so/fragments_2"
os.makedirs(dir_fragments, exist_ok=True)

# fragment BED 파일로 내보내기
#  - prefix="mcluster", suffix="_fragments.bed": 파일명 규칙
#  - 샘플×클러스터별로 파일이 생성됩니다.
snap.ex.export_fragments(
    data,
    groupby=groupby_series,   # ← 컬럼명 대신 시리즈를 직접 전달 (인덱스로 매핑)
    out_dir=dir_fragments,
    prefix="mcluster",
    suffix="_fragments.bed"
)


# ------------------------------------------------------------------------------------
# [9] *_updated.h5ad 파일을 다시 하나의 AnnData로 합치고 집계 통계 출력
#  - 경로: /aggr_update 폴더 (주의: 위에서 저장한 폴더명은 aggr_update_2 입니다)
#  - ad.concat으로 하나의 AnnData로 합쳐 크로스탭/카운트 출력
#  - 마지막에 하나짜리 h5ad(merged)로도 보관
# ------------------------------------------------------------------------------------
import os, glob
import scanpy as sc
import anndata as ad
import pandas as pd

# 실제 파일 있는 폴더로 설정 
base = "/home/sohyeong/projects/better_atac_pipeline/output_so/aggr_update"

# *_updated.h5ad 자동 수집
paths = sorted(glob.glob(os.path.join(base, "sample*_updated.h5ad")))
if not paths:
    raise FileNotFoundError(f"*_updated.h5ad 파일을 찾지 못했습니다: {base}")

# 읽고 합치기
ad_list, keys = [], []
for p in paths:
    a = sc.read_h5ad(p)
    # 필요한 obs 컬럼 확인/캐스팅
    for col in ["leiden", "celltype", "sample"]:
        if col not in a.obs:
            raise KeyError(f"{os.path.basename(p)}에 '{col}' 컬럼이 없습니다.")
        a.obs[col] = a.obs[col].astype(str)
    ad_list.append(a)
    keys.append(os.path.splitext(os.path.basename(p))[0].replace("_updated",""))

ad_all = ad.concat(ad_list, join="outer", label="sample", keys=keys, index_unique=None)

# 집계
obs = ad_all.obs
print("📊 클러스터별 cell 개수")
cluster_counts = obs["leiden"].value_counts().sort_index()
print(cluster_counts, "\n")

print("📊 세포타입별 cell 개수")
celltype_counts = obs["celltype"].value_counts()
print(celltype_counts, "\n")

print("📊 클러스터 × 세포타입 교차표")
cluster_celltype_table = pd.crosstab(obs["leiden"], obs["celltype"])
print(cluster_celltype_table, "\n")

print("📊 샘플 × 세포타입 교차표")
sample_celltype_table = pd.crosstab(obs["sample"], obs["celltype"])
print(sample_celltype_table)

# 📂 CSV로 저장
outdir = "/home/sohyeong/projects/better_atac_pipeline/output_so/aggr_update"
os.makedirs(outdir, exist_ok=True)

cluster_counts.to_csv(os.path.join(outdir, "cluster_counts.csv"))
celltype_counts.to_csv(os.path.join(outdir, "celltype_counts.csv"))
cluster_celltype_table.to_csv(os.path.join(outdir, "cluster_x_celltype.csv"))
sample_celltype_table.to_csv(os.path.join(outdir, "sample_x_celltype.csv"))

print(f"\n✅ CSV 파일로 저장 완료 (경로: {outdir})")