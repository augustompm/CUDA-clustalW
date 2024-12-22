#!/bin/bash

# Diretórios de destino
CLUSTALW_DATA_DIR="/dados/home/tesla-dados/CUDA-clustalW/data"
DEFAULT_SCORE_DIR="/dados/home/tesla-dados/CUDA-clustalW/default-score"

# Criação de diretórios, caso não existam
mkdir -p "$CLUSTALW_DATA_DIR"
mkdir -p "$DEFAULT_SCORE_DIR"

# URLs para download
BAliBASE_R1_5="https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz"
BAliBASE_R10="https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R10.tar.gz"
BAliBASE_R9_SCORE="https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R9_bali_score.tar.gz"
BAliBASE_R10_SCORE="https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R10_bali_score.tar.gz"

# Download dos arquivos
echo "Baixando arquivos..."
wget -q -P "$CLUSTALW_DATA_DIR" "$BAliBASE_R1_5"
wget -q -P "$CLUSTALW_DATA_DIR" "$BAliBASE_R10"
wget -q -P "$DEFAULT_SCORE_DIR" "$BAliBASE_R9_SCORE"
wget -q -P "$DEFAULT_SCORE_DIR" "$BAliBASE_R10_SCORE"

# Extração dos arquivos
echo "Extraindo arquivos..."
tar -xzf "$CLUSTALW_DATA_DIR/BAliBASE_R1-5.tar.gz" -C "$CLUSTALW_DATA_DIR"
tar -xzf "$CLUSTALW_DATA_DIR/BAliBASE_R10.tar.gz" -C "$CLUSTALW_DATA_DIR"
tar -xzf "$DEFAULT_SCORE_DIR/BAliBASE_R9_bali_score.tar.gz" -C "$DEFAULT_SCORE_DIR"
tar -xzf "$DEFAULT_SCORE_DIR/BAliBASE_R10_bali_score.tar.gz" -C "$DEFAULT_SCORE_DIR"

# Remoção dos arquivos .tar.gz
echo "Limpando arquivos .tar.gz..."
rm "$CLUSTALW_DATA_DIR/BAliBASE_R1-5.tar.gz"
rm "$CLUSTALW_DATA_DIR/BAliBASE_R10.tar.gz"
rm "$DEFAULT_SCORE_DIR/BAliBASE_R9_bali_score.tar.gz"
rm "$DEFAULT_SCORE_DIR/BAliBASE_R10_bali_score.tar.gz"

echo "Arquivos organizados com sucesso."

