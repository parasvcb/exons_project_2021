pip install Django==4.1.2
pip install djangorestframework
pip install whitenoise


mysqldump -u root -p enactdb > enactd.sql

MariaDB [(none)]> show databases;
+--------------------+
| Database           |
+--------------------+
| information_schema |
| mysql              |
| performance_schema |
+--------------------+
3 rows in set (0.045 sec)

MariaDB [(none)]> CREATE DATABASE enactdb;
Query OK, 1 row affected (0.000 sec)

MariaDB [(none)]> Exit;
Bye


# installing mysql has been explained in the component 2 

after changing the models and serilaizers, 
do 

python manage.py migrate --fake <appname>
and then
python manage.py runserver



'''
# script that builds reactapp and moves it to our zone
#~/bin/bash
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
djangoDir=$(cd $djangoDir; pwd)
rm -rf ${djangoDir}/frontend #$djangoDir=/var/www/html/nextrap2/
mkdir -p ${djangoDir}/frontend/build/root
rsync -a build/ ${djangoDir}/frontend/build/ && cd ${djangoDir}/frontend/build
mv *.ico *.js *.json root/
cd $djangoDir
rm -rf static && mkdir static
python manage.py collectstatic
'''
