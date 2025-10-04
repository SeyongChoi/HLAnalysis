 HLAnalysis

**HLAnalysis**는 Solid/Water interface 시스템에서 solid-water interface에서 형성되는 **Hydration Layer(수화층)**의 LAMMPS trajectory로부터 원자 구조를 읽고 이를 통해 구조·동역학·스펙트럼 분석을 구현한 프로젝트입니다. 이 프로젝트는 Rust 및 python 언어를 기반으로 하며, density profile, molecular orientation, time orientational correlation function, SFG spectrum을 포함하여 다양한 실험 구성이 가능하도록 설계되어 있습니다.


## 📁 프로젝트 구조

```
HLAnalysis/
├── pyproject.toml              # maturin 기반 Python 빌드 설정 (PEP 621)
├── Cargo.toml                  # Rust 크레이트 설정
├── src/                        # Rust 소스코드 (핵심 로직)
│   ├── lib.rs                  # PyO3 엔트리 포인트
│   ├── atoms/                  # 원자 데이터 구조
│   │   ├── mod.rs
│   │   ├── atom.rs
│   │   └── atoms.rs
│   ├── io/                     # 입력/출력 (I/O)
│   │   ├── mod.rs
│   │   ├── read/
│   │   │   ├── mod.rs
│   │   │   ├── read.rs
│   │   │   ├── read_xyz.rs
│   │   │   └── read_lammps_dump.rs
│   │   └── write/
│   │       ├── mod.rs
│   │       ├── write.rs
│   │       ├── write_hbond_info.rs
│   │       └── write_den_profiles.rs
│   ├── analysis/               # 분석 루틴 (핵심 계산 모듈)
│   │   ├── mod.rs
│   │   ├── density_profile/
│   │   │   ├── mod.rs
│   │   │   ├── density_profile.rs
│   │   │   └── README.md
│   │   ├── orient_dist/
│   │   │   ├── mod.rs
│   │   │   ├── bond_orient_along_normal.rs
│   │   │   ├── molecular_orient_dist.rs
│   │   │   └── README.md
│   │   ├── orient_dynamics/
│   │   │   ├── mod.rs
│   │   │   ├── orient_time_corr_function.rs
│   │   │   └── README.md
│   │   └── sfg/
│   │       ├── mod.rs
│   │       ├── dipole_polar_term.rs
│   │       ├── vvacf.rs
│   │       ├── tr_vvacf.rs
│   │       ├── spectrum.rs
│   │       └── README.md
│   └── utils/                  # 보조 유틸리티 모듈
│       ├── mod.rs
│       ├── parameters.rs
│       ├── atomic_data.rs
│       ├── mic.rs
│       ├── converter.rs
│       └── pycell3.rs
├── python/                     # Python wrapper 및 후처리 모듈
│   └── hlanalysis/
│       ├── __init__.py         # Rust 모듈 래핑 및 네임스페이스 초기화
│       ├── _version.py         # 버전 정보
│       ├── io.py               # Rust 함수 래퍼 또는 파일 입출력 후처리
│       ├── analysis/           # 후처리 및 시각화 (Rust 결과를 Python에서 가공)
│       │   ├── __init__.py
│       │   ├── plot_density.py
│       │   ├── plot_orient.py
│       │   └── plot_sfg.py
│       └── utils/
│           ├── __init__.py
│           ├── filepaths.py
│           └── timer.py
├── requirements.txt             # Python 의존성 목록
└── README.md                    # 프로젝트 설명 문서

```
<!-- 
## 🧠 주요 기능

- LAMMPS trajectory로 부터 atomic information 읽기
- Density profile
- Anaylsis molecular orientation
- Orientational time correlation function
- SFG spectrum
- 결과 시각화 자동화

## ⚙️ 설치 방법

```bash
# Conda 환경 예시
conda create -n steerablecnnca python=3.8
conda activate steerablecnnca

# 저장소 클론
git clone https://github.com/SeyongChoi/SteerableCNNCA.git
cd SteerableCNNCA

# 필수 패키지 설치
pip install -r requirements.txt
```

## 🚀 실행 방법

### 1. 설정 파일 준비

`config/` 디렉토리의 YAML 파일을 수정하여 데이터 경로, 모델 종류, 학습 설정 등을 구성합니다.

```yaml
model:
  type: "SteerableCNN"  # or "CNN", "ANN"

dataset:
  data_root_dir: "./data/"
  grid_size: 100
  ...
```

### 2. 학습 실행

```bash
python main.py --config config/steerablecnn.yaml
```

학습 로그는 W&B와 콘솔에 출력되며, 모델 체크포인트 및 예측 결과는 지정된 출력 폴더에 저장됩니다.

## 🧩 모델 종류

- `ANNModel`: MLP 기반 단순 회귀
- `CNNModel`: 2D ConvNet 기반 회귀
- `SteerableCNNModel`: 회전 변환에 불변한 steerable filter 기반 CNN (ESCNN 사용)

## 📊 예측 결과 예시

- 학습 loss curve  
- 예측값 vs 실제값 scatter plot  
- 격자 데이터 시각화

<p align="center">
  <img src="docs/example_plot.png" width="500">
</p>

## 📦 주요 의존성

- Python 3.8+
- PyTorch
- PyTorch Lightning
- [ESCNN](https://github.com/QUVA-Lab/escnn)
- wandb
- numpy, matplotlib, scikit-learn 등

## ✍️ 작성자

- **Seyong Choi** – [GitHub 프로필](https://github.com/SeyongChoi)

## 📄 라이선스

본 프로젝트는 MIT 라이선스를 따릅니다. (필요 시 명시)

---

이 문서는 학습, 실험, 평가를 효율적으로 관리하고자 하는 사용자와 연구자를 위한 안내서입니다.  
피드백이나 제안이 있다면 언제든지 Issue나 PR을 통해 공유해주세요!
 -->