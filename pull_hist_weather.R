
startdate = as.Date("2014/01/01")
enddate = as.Date("2014/06/30")

firstpart = "http://www.wunderground.com/history/airport/EGFF/"
lastpart = "/DailyHistory.html?req_city=Cardiff&req_state=&req_statename=United+Kingdom&reqdb.zip=00000&reqdb.magic=1&reqdb.wmo=03717&format=1"

##weather_data = read.csv(paste(firstpart,date,lastpart)) 
hist_weather_data = list()
for (i in 0:(enddate-startdate)) {  
  
  date2 = gsub("-","/", startdate+i)
  hist_weather_data[[i+1]] = read.csv(paste(firstpart,date2,lastpart,sep="")) 
  
}

hist_wind = unlist(sapply(hist_weather_data,function(x){x$Wind.SpeedMPH},simplify = TRUE))
hist_pres = unlist(sapply(hist_weather_data,function(x){x$Sea.Level.PressureIn},simplify = TRUE))
hist_pres[hist_pres < 5] = NA
dateUTC = unlist(sapply(hist_weather_data,function(x){x$DateUTC},simplify = TRUE))

day = strptime(sub(pattern = "<br />","",as.character(dateUTC)),format="%Y-%m-%d %T")

#met = Reduce(function(...) merge(...,all=TRUE),weather_data)

##jpeg(file="~/SURGE/web/WindSpeed.jpg")
plot(day,hist_wind, ylab="Wind Speed (mph)", xlab="Date", main="Historical Wind Speed in 30 min Intervals in Cardiff,UK",type='l')
##dev.off()

##jpeg(file="~/SURGE/web/SeaLevelPressure.jpg")
plot(day,hist_pres, ylab="Sea Level Pressure (in)", xlab="Date", main="Historical Sea Level Pressure in 30 min Intervals in Cardiff,UK",type='l')
##dev.off()