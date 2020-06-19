#!/bin/bash -l

# Exit immediately on errors
set -e

function echo_run {
	echo "$" "$@"
	"$@"
}

REPO_SRC="${GITHUB_WORKSPACE}/${INPUT_SOURCE_DIR}"
GH_PAGES="${GITHUB_WORKSPACE}/${INPUT_PAGES_DIR}"
SPHNX_DT=$GITHUB_WORKSPACE/.doctree

echo ::group::Create working directories
echo_run mkdir -p $SPHNX_DT
HTML_DIR=`mktemp -d -p ${GITHUB_WORKSPACE}`
echo ::endgroup::

echo ::group::Gather author and commit information
echo_run cd $REPO_SRC
AUTHOR_NAME="$(git show --format=%an -s)"
AUTHOR_EMAIL="$(git show --format=%ae -s)"
DOCS_SHA8="$(echo ${GITHUB_SHA} | cut -c 1-8)"
: "${named:=$(git describe --tags)}" "${named:=$(git describe --all)}"
echo "::set-output name=name::"${AUTHOR_NAME}""
echo "::set-output name=email::"$AUTHOR_EMAIL}""
echo "::set-output name=docs_sha::$(echo ${GITHUB_SHA})"
echo "::set-output name=docs_sha8::"${DOCS_SHA8}""
echo "git described as ${named}"
echo ::endgroup::

echo ::group::Build and install the python module
#echo_run python3 -m pip install $REPO_SRC
echo_run python3 -m pip install brille
echo ::endgroup::

echo ::group::Configure pages author information
echo_run cd $GH_PAGES
echo_run git config user.name $AUTHOR_NAME
echo_run git config user.email $AUTHOR_EMAIL
echo ::endgroup::

echo ::group::Build pages with Sphinx
# pre-create the _build directory so that Doxygen is happy (This is the opposite of good)
mkdir -p $REPO_SRC/docs/_build
# actually build the HTML pages
echo_run sphinx-build -b html $REPO_SRC/docs $HTML_DIR -E -d $SPHNX_DT
echo ::endgroup::

if [ "${INPUT_CREATE_README}" = true ]; then
	echo ::group::Create README.md
	README_FILE="${HTML_DIR}/README.md"
	URL_BASE="https://github.com/${GITHUB_REPOSITORY}"
	echo_run touch $README_FILE
	echo "GitHub Pages of [${GITHUB_REPOSITORY}](${URL_BASE}.git)" > $README_FILE
	echo "======================================" >> $README_FILE
	echo "Sphinx HTML documentation of [${DOCS_SHA8}](${URL_BASE}/tree/${GITHUB_SHA})" >> $README_FILE
	echo ::endgroup::
fi


echo ::group::Move documentation into the correct places
named_dir="${GH_PAGES}/latest"
if [ "${INPUT_IS_RELEASE}" = true ]; then
	# Cleanup the gh-pages root
	echo_run cd $GH_PAGES
	echo_run git rm -rf --quiet .
	# unsage delection and recover directories
	echo "git reset then checkout stable and latest directories"
	git reset --quiet -- stable || true
	git reset --quiet -- latest || true
	git checkout --quiet -- stable || true
	git checkout --quiet -- latest || true
	# Copy the newly built HTML pages into the gh_pages root
	echo_run rsync -a "${HTML_DIR}/" $GH_PAGES
	echo_run touch "${GH_PAGES}/.nojekyll"
	# Switch from gh_pages/latest to (ideally) gh_pages/stable/[version name]
	named_dir="${GH_PAGES}/stable/${named}"	
fi
# Copy the newly built HTML pages into either gh_pages/latest or gh_pages/stable/[version name]
echo_run mkdir -p $named_dir
echo_run rsync -a --delete "${HTML_DIR}/" $named_dir
# Add all changes to the gh_pages branch
echo_run cd $GH_PAGES
echo_run git add .
echo ::endgroup::

echo ::group::Commit to repository, overwriting last commit, and push
echo_run git commit --amend --date="$(date)" --message='Auto commit from Versioned Sphinx GH-Pages Action'
echo_run git push --force-with-lease
echo ::endgroup::
