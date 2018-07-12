
#First install the reporter tool
wget -O ~/codacy-coverage-reporter-assembly-latest.jar https://oss.sonatype.org/service/local/repositories/releases/content/com/codacy/codacy-coverage-reporter/4.0.0/codacy-coverage-reporter-4.0.0-assembly.jar

#Analyse the existing gcov files and output in compatible xml format
#The -g option says to use the existing gcov files, the -k options says to not delete the gcov files
#The -j option is like make's -j
gcovr --root . -g -k -j 2 --xml -o gcovr_report.xml

#Do the upload
java -jar ~/codacy-coverage-reporter-assembly-latest.jar report --language CPP -r gcovr_report.xml --partial
