#!/bin/bash -l

# Exit immediately on errors
set -e

function echo_run {
	echo "$" "$@"
	"$@"
}

REPO_SRC="${GITHUB_WORKSPACE}/${INPUT_SOURCE_DIR}"
SPHNX_DT=$GITHUB_WORKSPACE/.doctree

echo ::group::Create working directories
echo_run mkdir -p $SPHNX_DT
# Create a temporary directory *in* the Docker container to avoid cleaning-up problems
HTML_DIR=`mktemp -d`
echo ::endgroup::

echo ::group::Gather author and commit information
echo_run cd $REPO_SRC
AUTHOR_NAME="$(git show --format=%an -s)"
AUTHOR_EMAIL="$(git show --format=%ae -s)"
DOCS_SHA8="$(echo ${GITHUB_SHA} | cut -c 1-8)"

BRANCH=$(git branch --show-current)
if [ $BRANCH = "master" ]; then
  # Determine the tag name or branch name -- nevermind, tag or 'latest'
  # only match 'version' tags, e.g., vM.m.p, and keep only up to the second period
  : "${named:=$(git describe --tags --match 'v*'| cut -d'.' -f1-2)}" "${named:=latest}"
else
  named="branch-$BRANCH"
fi

echo "::set-output name=name::"${AUTHOR_NAME}""
echo "::set-output name=email::"$AUTHOR_EMAIL}""
echo "::set-output name=docs_sha::$(echo ${GITHUB_SHA})"
echo "::set-output name=docs_sha8::"${DOCS_SHA8}""
echo "git described as ${named}"
echo ::endgroup::

echo ::group::Build and install the python module
# echo_run python3 -m pip install --use-feature=in-tree-build $REPO_SRC
echo_run cd $REPO_SRC
echo_run python3 setup.py install
echo ::endgroup::


echo ::group::Build pages with Sphinx
# create directories to make Doxygen happy (This is the opposite of good)
mkdir -p $REPO_SRC/docs/_build/doxygenxml
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

for val in $INPUT_PAGES_DIR; do
	echo ::group::Move documentation into the correct places
	named_dir="${GITHUB_WORKSPACE}/${val}/${named}"
	# Copy the newly built HTML pages into either gh_pages/latest or gh_pages/[version name]
	echo_run mkdir -p $named_dir
	echo_run rsync -a --delete "${HTML_DIR}/" $named_dir
	# make sure we're in the root of the gh-pages branch
	echo_run cd "${GITHUB_WORKSPACE}/${val}"
	# for releases, move the stable symlink
	if [ "${INPUT_IS_RELEASE}" = true ]; then
		echo_run unlink stable
		echo_run ln -s ${named} stable
	fi
	echo ::endgroup::
  if [ "${INPUT_UPDATE_GIT}" = true ]; then
    echo ::group::Configure pages author information
    echo_run git config user.name $AUTHOR_NAME
    echo_run git config user.email $AUTHOR_EMAIL
    echo ::endgroup::
	  echo ::group::Add all changes to the repository
	  echo_run git add .
	  echo ::endgroup::
	  echo ::group::Commit to repository and push
	  echo_run git commit --date="$(date)" --message='Auto commit from Versioned Sphinx GH-Pages Action'
	  echo_run git push
	  echo ::endgroup::
  fi
done
