bundle exec jekyll clean --force
rm -rf _site
bundle clean --force
bundle install
bundle exec jekyll build
