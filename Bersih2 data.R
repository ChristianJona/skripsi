library(tidyverse)

# Path ke file CSV
path_females <- "C:/Users/rapha/Documents/1 PROJECT/0 SKRIPSI/Raw Data dan R/Data/Female Belanda.csv"
path_males   <- "C:/Users/rapha/Documents/1 PROJECT/0 SKRIPSI/Raw Data dan R/Data/Male Belanda.csv"

# Baca data FEMALES
data_females <- read_csv2(path_females) %>%
  mutate(Sex = "Female") %>%
  rename(Cohort = Year)

# Baca data MALES
data_males <- read_csv2(path_males) %>%
  mutate(Sex = "Male") %>%
  rename(Cohort = Year)

# Gabungkan dan bersihkan
data_gabungan <- bind_rows(data_females, data_males)

data_dibersihkan <- data_gabungan %>%
  mutate(
    Age = as.numeric(Age),
    lx  = as.numeric(lx),
    dx  = as.numeric(dx),
    qx  = as.numeric(qx)
  ) %>%
  filter(!is.na(Age))

# Filter nilai relevan
data_siap_analisis <- data_dibersihkan %>%
  filter(
    Cohort >= 1893,
    Cohort <= 1908,
    Age >= 65
  ) %>%
  select(Cohort, Age, Sex, lx, dx, qx)

# Cek hasil
cat("Data Siap untuk Dianalisis (setelah paksa konversi tipe data):\n")
cat("Jumlah baris:", nrow(data_siap_analisis), "\n")
cat("Rentang umur:", range(data_siap_analisis$Age, na.rm = TRUE), "\n")
cat("Rentang kohor:", range(data_siap_analisis$Cohort, na.rm = TRUE), "\n\n")

print(head(data_siap_analisis))
print(tail(data_siap_analisis))

str(data_siap_analisis)

# Simpan sebagai RDS
saveRDS(data_siap_analisis,
        file = "C:/Users/rapha/Documents/1 PROJECT/0 SKRIPSI/Raw Data dan R/data bersih/data_final_siap_pakai.rds")
