set -ex

if which docker &> /dev/null ; then
    cmd="time sudo docker"
else
    cmd="time podman"
fi

file=$1
(test $file && test -f $file) || file=.docker/fedora/Dockerfile
test $# -gt 0 && shift

if test -x $file
then
    COMMIT=$(git rev-parse HEAD) $file $@ > Dockerfile
else
    cp $file Dockerfile
fi

$cmd login -p $DOCKER_TOKEN -u $DOCKER_USER

cat Dockerfile

$cmd build -t mobydick . --build-arg=COMMIT=$(git rev-parse HEAD) --build-arg=MPI=$MPI --build-arg=CMAKE_OPTIONS="$CMAKE_OPTIONS"
for tag in $TAGS
do
    $cmd tag mobydick $tag
    $cmd push $tag
done
