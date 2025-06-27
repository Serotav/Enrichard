SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MODULES_DIR="$SCRIPT_DIR/../app/Modules"
WEHI_TEMPLATE="$SCRIPT_DIR/wehi_template"
COMMON_BACKGROUND="$SCRIPT_DIR/../app/Common_Background"

wehi_sets=(
#"Hallmark_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds" 
"Positional_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c1.all.v7.1.entrez.rds" 
"Curated_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.all.v7.1.entrez.rds"
#"Regulatory_Target_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c3.all.v7.1.entrez.rds"
#"Computational_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c4.all.v7.1.entrez.rds"
#"Ontology_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c5.all.v7.1.entrez.rds"
#"Oncogenic_Signature,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c6.all.v7.1.entrez.rds"
#"Immunologic_Signature_Gene,https://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c7.all.v7.1.entrez.rds"
)

echo $COMMON_BACKGROUND
rm "$COMMON_BACKGROUND"/*

rm -rf "$MODULES_DIR"
mkdir -p "$MODULES_DIR"

for item in "${wehi_sets[@]}"; do
    IFS=',' read -r set url <<< "$item"
    printf "%-30s %s\n" "$set" "$url"
    cp -r $WEHI_TEMPLATE "$MODULES_DIR/$set"
    echo "$url" > "$MODULES_DIR/$set/source_url.txt"
done

