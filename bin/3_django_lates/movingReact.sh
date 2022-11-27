#~/bin/bash
# real file is in ~/bin or production server
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, please enter the 1. reactDir 2. staticDir (db.sqlite3 exonapp exon.db frontend manage.py static v1)"
    exit 0
fi

reactDir=$1
djangoDir=$2
djangoDir=$(cd $djangoDir; pwd)
reactDir=$(cd $reactDir; pwd)
#pwd
#echo $djangoDir
#echo $reactDir
#exit 0
set python="~/miniconda3/envs/django4/bin/python"
cd $reactDir
npm run build
grep -rl "http://localhost:3000" build/ | xargs sed -i -e 's,http://localhost:3000,http://172.16.46.73/enactdb,g'
grep -rl "http://localhost:8000" build/ | xargs sed -i -e 's,http://localhost:8000,http://172.16.46.73/enactdb,g'
djangoDir=$(cd $djangoDir; pwd)
rm -rf ${djangoDir}/frontend #$djangoDir=/var/www/html/nextrap2/
mkdir -p ${djangoDir}/frontend/build/root
rsync -a build/ ${djangoDir}/frontend/build/ && cd ${djangoDir}/frontend/build
mv *.ico *.js *.json root/
cd $djangoDir
rm -rf static && mkdir static
python manage.py collectstatic

~                                                                                                                                                 
~                                                                                                                                                 
~                                                                                                                                             
