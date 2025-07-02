# Variant Frequency Calculation Workflow

Bu klasör, modern ve karşılaştırmalı örneklerde varyant çağrımı yapıp, fenotip/grup bilgileri ekleyerek her grup için risk/minor alel frekanslarını maksimum olasılık yöntemini kullanarak hesaplayan **üç aşamalı** bir iş akışını içerir.

> **Akış Özeti**
>
> 1. **`variant_call.sh`** » Örnek BAM dosyalarından VCF oluşturur.
> 2. **`process_variant_table.R`** » Varyant tablosuna fenotip, risk/minor alel ve grup bilgileri ekler.
> 3. **`freq_calculation.R`** » Her grup‑SNP‑fenotip kombinasyonu için alel frekansını ve %95 Wilson güven aralığını hesaplar.

---

## 1 | Ön Gereksinimler

| Yazılım                                                                                                                                               | Versiyon / Not |
| ----------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| **bcftools**                                                                                                                                          | ≥ 1.18         |
| Yol, `variant_call.sh` içinde `bcftools=` satırında tanımlı                                                                                           |                |
| **R**                                                                                                                                                 | ≥ 4.2          |
| İlave paketler: `tidyverse`, `data.table`, `ggpubr`, `dplyr`, `Hmisc`, `glue`, `binom`                                                                |                |
| **SLURM**                                                                                                                                             | (İsteğe bağlı) |
| `variant_call.sh` SLURM job direktifleri içerir; yerel bilgisayarda çalıştırılacaksa `#!/bin/bash` kısmı dışındaki `#SBATCH` satırları kaldırılabilir |                |

Ayrıca:

* **Referans genom** (hg38 FASTA)
* **Bölge BED dosyası** (50 kb aralıkları)
* **Örnek BAM dizini**
* **Örnek isimleri listesi** (`alkan_samples.txt`, her satırda bir BAM adı)
* **Risk\_allele** tablosu
  (Sütunlar: *Chromosome, POS, Risk\_Allele, Minor\_Allele, …, Phenotype*)

---

## 2 | Dizin Yapısı

```
frequency_calculation/
├── variant_call.sh
├── process_variant_table.R
├── freq_calculation.R
└── vcf_files/          # Çıktı VCF'ler buraya düşer (otomatik oluşturulur)
```

---

## 3 | Adım Adım Çalıştırma

### 3.1  Varyant Çağrımı — `variant_call.sh`

```bash
# SLURM’da göndermek için
sbatch variant_call.sh

# Yerel çalıştırmak için (SBATCH satırlarını silin)
bash variant_call.sh
```

Script neler yapar?

1. **Değişkenler**: `bcftools`, `ref` (hg38 FASTA), `bedfile`, `INFILE`, `OUTDIR` tanımlanır.
2. `OUTDIR` yoksa oluşturulur → `vcf_files/`.
3. `INFILE` listesindeki her örnek için:

   * `bcftools mpileup` → DP & AD field’ları dâhil, BED bölgesiyle sınırlı mpileup akışı.
   * `bcftools call -mv` → VCF üretilir (`${sample}.vcf`).
4. Her VCF’e örnek adı eklenir (`awk ... > .vcf.modified`).
5. Tüm modifiye VCF’ler **`combined.vcf`**’te birleştirilir.

> **Çıktı**: `vcf_files/*.vcf` + `combined.vcf`

### 3.2  Varyant Tablosu İşleme — `process_variant_table.R`

```R
# Etkileşimli R oturumunda veya
Rscript process_variant_table.R
```

Script aşağıdakileri yapar:

1. Sunucudan gelen CSV içe alınır (sütun başlıkları: `CHROM, POS, DP, DP4, REF, ALT, BamID`).
2. Pozisyon bilgilerinden `POS` kimliği oluşturur (`chr_pos`).
3. `Risk_allele` tablosuyla eşleştirip:

   * `SNP_ID`, `Risk_Allele`, `Minor_Allele`, `Phenotype` ekler.
4. `DP4` alanını dört kolona ayırır (`Ref1,Ref2,Alt1,Alt2`).
5. **Risk (veya Minor) alel ok sayısı** hesaplanır ➜ `Total_RAC`.
6. **Grup** bilgisi eklenir (örnekte tümü `"Modern"`; tarihsel vs. gruplar için düzenleyin).
7. Nihai çıktı veri çerçevesi → `modern_df` (sütunlar: `SampleName, Group, POS, SNP_ID, Phenotype, R, T`).

> **Not**: Dosya yolları ve nesne adları (ör. `comperative_alkan`, `Risk_allele`) betikte sabit; kendi veri adlarınıza göre güncelleyin.

### 3.3  Frekans Hesaplama — `freq_calculation.R`

```R
Rscript freq_calculation.R
```

1. `modern_df` ve (opsiyonel) başka grup veri çerçeveleri (`comp_others4` vb.) birleştirilir.
2. Her **Grup × SNP\_ID × Phenotype** alt‑kümesi için:

   * Maksimum Olasılık Tahmini (MLE) ile risk/minor alel frekansı `p̂` bulunur.
   * Wilson %95 güven aralığı hesaplanır (`binom.confint`).
3. Sonuçlar → **`all_freq.csv`** (varsayılan yol: masaüstü). Sütunlar:

   * `Period (Group), SNPid, Phenotype, POPsize, pHat, CI_Lower, CI_Upper`, vs.

> **Çıktı**: `all_freq.csv`

---

## 4 | Sonuç Dosyaları

| Dosya                  | Açıklama                                        |
| ---------------------- | ----------------------------------------------- |
| `combined.vcf`         | Tüm örneklerin tek VCF dosyası                  |
| `modern_df` (R objesi) | Fenotip ve grup eklenmiş varyant tablosu        |
| `all_freq.csv`         | Grup bazlı alel frekansları ve güven aralıkları |

---

## 5 | Özelleştirme & İpuçları

* **SLURM direktifleri**: Küme kaynaklarınıza göre `--partition`, `--time`, `--nodes` vb. alanları güncelleyin.
* **Bölge seçimi**: `bedfile` satırını değiştirerek farklı genom bölgeleri hedefleyin.
* **Grup Ataması**: `process_variant_table.R` içinde `modern$Group = "Modern"` satırını kendi sınıflandırmanıza göre düzenleyin.
* **Rscript Parametreleştirme**: Betikler şu anda sabit yol/nesne adlarıyla yazılmıştır; argümanla değişken almak için `commandArgs(trailingOnly=TRUE)` eklenebilir.

---

## 6 | Referanslar

* **bcftools**: Danecek et al., 2011, *Bioinformatics* 27(21): 2987‑2993.
* **Wilson Güven Aralığı**: Wilson, 1927, *J. Amer. Stat. Assoc.* 22: 209‑212.

---

> Soru, yorum veya katkılarınız için lütfen *Issues* bölümünden çekinmeden ulaştırın.
