set -ex

if which docker &> /dev/null ; then
    cmd="time sudo docker"
else
    cmd="time podman"
fi

file=$1
(test $file && test -f $file) || file=.docker/fedora/Dockerfile.sh
test $# -gt 0 && shift
if test -e $file
then
    COMMIT=$(git rev-parse HEAD) bash $file $@ > Dockerfile
else
    cp $file Dockerfile
fi

$cmd login -p $DOCKER_TOKEN -u $DOCKER_USER

cat Dockerfile

$cmd build -t mobydick .
for tag in $TAGS
do
    $cmd tag mobydick $tag
    $cmd push $tag
done
