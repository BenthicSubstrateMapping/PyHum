bundle exec jekyll clean
rm -rf _site
bundle clean
bundle install
bundle exec jekyll build
