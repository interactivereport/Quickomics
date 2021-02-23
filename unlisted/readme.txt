This is the folder where Rdata generated from user uploaded csv files will be saved.
Please don't delete the unlisted/ folder, otherwise the Upload Files tool will not work.

The Upload Files tools from QuickOmics system will create three files for each datasets in the unlisted/ folder
*ProjectID.RData  (The main R data files)
*ProjectID_network.RData  (gene-gene correlation network data)
*ProjectID.csv (the same format as the saved_projects.csv in data/ folder) 


One can also upload pre-made RData files to this folder, and allow users to access the projet using unlisted URL: 
http://127.0.0.1:5772/?unlisted=ProjectID 
Replace 127.0.0.1:5772 with the actualy URL of the R Shiny App. Replace ProjectID with the actual ID of the unlisted project.
Unlist projects are only visible to users with the special URL. 
