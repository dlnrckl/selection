#!/bin/bash

# Dosya adı, satır sayısı ve .bed dosyası oluşturmak için
while read -r line; do
    # Satırdaki kolonları ayır
    set -- $line
    # İlgili verileri al
    chrom=$1
    start=$2
    end=$3
    name=$4
    
    # .bed uzantılı dosyaya yaz
    echo -e "$chrom\t$start\t$end\t$name" > "${name}.bed"
done < 50kb_hg38
