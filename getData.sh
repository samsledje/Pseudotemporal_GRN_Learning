mkdir -p data
curl -L -o tmpData.zip https://www.dropbox.com/sh/anczobj20kiyyqf/AAANJhHO8pL_hVM329EeHMkTa?dl=1
mv tmpData.zip data
cd data
unzip tmpData.zip
rm tmpData.zip
cd ..
