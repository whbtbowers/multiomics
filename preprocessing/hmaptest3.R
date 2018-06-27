library(ggplot2)
#Reading in the data
chicagoMVT <- read.csv('../data/motor_vehicle_theft.csv', stringsAsFactors = FALSE)
#Converting the date to a recognizable format
chicagoMVT$Date <- strptime(chicagoMVT$Date, format = '%m/%d/%Y %I:%M:%S %p')
#Getting the day and hour of each crime
chicagoMVT$Day <- weekdays(chicagoMVT$Date)
chicagoMVT$Hour <- chicagoMVT$Date$hour
#Sorting the weekdays
dailyCrimes <- as.data.frame(table(chicagoMVT$Day, chicagoMVT$Hour))
names(dailyCrimes) <- c('Day', 'Hour', 'Freq')
dailyCrimes$Hour <- as.numeric(as.character(dailyCrimes$Hour))
dailyCrimes$Day <- factor(dailyCrimes$Day, ordered = TRUE, 
                          levels = c('Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'))
#Plotting the number of crimes each day (line graph)
ggplot(dailyCrimes, aes(x = Hour, y = Freq)) + geom_line(aes(group = Day, color = Day)) + xlab('Hour') + ylab('Number of thefts') + ggtitle('Daily number of Motor Vehicle Thefts')
