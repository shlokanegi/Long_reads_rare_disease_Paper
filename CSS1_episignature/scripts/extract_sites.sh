
tabix $1 -R data/CSS1.cpgs.p01.regions > $(echo $1 | sed 's/beds/site_beds/'  | sed 's/.gz//')
