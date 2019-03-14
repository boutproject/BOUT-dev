#First install the reporter tool
wget -q -O ~/codacy-coverage-reporter-assembly-latest.jar https://oss.sonatype.org/service/local/repositories/releases/content/com/codacy/codacy-coverage-reporter/4.0.5/codacy-coverage-reporter-4.0.5-assembly.jar > /dev/null

# Create all the .gcov files
# Run gcov in the same directory as the source file for the main
# library (because we used recursive make)
find . -type f -name *.gcno -execdir gcov -pb -r {} +
# but just in the unit tests top-level dir as we ran make from there
find tests/unit -type f -name *.gcno -exec gcov -pb -r {} +

#Analyse the existing gcov files and output in compatible xml format
#The -g option says to use the existing gcov files, the -k options says to not delete the gcov files
#The -j option is like make's -j
#See https://gcovr.com/guide.html for more details
gcovr --root . -k -j 2 --xml -o gcovr_report.xml --exclude-directories "tests/.*"

#Do the upload
java -jar ~/codacy-coverage-reporter-assembly-latest.jar report -f --language CPP -r gcovr_report.xml --partial
