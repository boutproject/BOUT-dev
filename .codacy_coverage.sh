
#First install the reporter tool
wget -O ~/codacy-coverage-reporter-assembly-latest.jar $(curl https://api.github.com/repos/codacy/codacy-coverage-reporter/releases/latest | jq -r .assets[0].browser_download_url)

#Analyse the existing gcov files and output in compatible xml format
#The -g option says to use the existing gcov files, the -k options says to not delete the gcov files
#The -j option is like make's -j
gcovr --root . -g -k -j 2 --xml -o gcovr_report.xml

#Do the upload
java -jar ~/codacy-coverage-reporter-assembly-latest.jar report --language CPP --forceLanguage -r gcovr_report.xml --partial
