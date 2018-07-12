TODAY=$(date +"%Y-%m-%d")
mkdir -p ../img/$TODAY

read -n1 -p  "Extract patients by ID? (y/n) " EXTRACT

#echo $EXTRACT

if [ $EXTRACT == "n" ];
then
    echo -e "\n\nExisting extraction files will be used."
elif [ $EXTRACT == "y" ];
then
    echo -e "\n\nCreating new extraction files"
fi
