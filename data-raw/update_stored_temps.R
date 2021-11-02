# Description: ------------------------------------------------------------
#
# Update saved station data and write to extdata folder
#  - run every 3-6 months to keep number of requested rows in cliflo account low
#
# Built under R version 4.1.0 (2021-05-18)
# Simon Howard; howards@landcareresearch.co.nz | si.w.howard@gmail.com
#-------------------------------------------------------------------------#

# get user name and password from system environment
username <- Sys.getenv("cliflo_usrid")
password <- Sys.getenv("cliflo_pwd")


## CliFro request latest temperatures

# check CliFlo credentials
if(is.null(username) | is.null(username)) stop("Include CliFlo login and password when calling append_TempSeq(); see 'https://cliflo.niwa.co.nz/'")

# add user & stations
me <- clifro::cf_user(username = username,
                      password = password)
my.dts <- clifro::cf_datatype(select_1 = 4, select_2 = 2, check_box = 1)
# 15752 - Musselburgh station ID - # 5893 - Nugget Point Aws
my.stations = clifro::cf_station(15752, 5893)

# date accessed
access.date <- Sys.Date()

# fetch data
cf.temp_latest <- clifro::cf_query(user = me,
                                   datatype = my.dts,
                                   station = my.stations,
                                   start_date = strftime("2000-01-01", format = "%Y-%m-%d 00"),
                                   end_date = strftime(Sys.Date(), format = "%Y-%m-%d 00"))

# set access.date attribute
attr(cf.temp_latest, "access.date") <- Sys.Date()

# write to csv extdata folder package install directory
write.csv(cf.temp_latest, paste0("inst/extdata/saved_station_temp_", format(access.date, "%Y%m%d"), ".csv"), row.names = FALSE)

saved_station_temps <- cf.temp_latest

# also save as internal dataset
usethis::use_data(saved_station_temps)
