motor = read.csv("Documents/real data/clean/motor.csv", header=T, fileEncoding="big-5")
car = read.csv("Documents/real data/clean/car.csv", header=T, fileEncoding="big-5")
head(car)

con = c("新北市永和區", "新北市蘆洲區", "臺北市大安區", "高雄市新興區", "新北市板橋區",
        "新北市三重區", "新北市新莊區", "臺北市大同區", "臺北市松山區", "臺中市北區",
        "高雄市苓雅區", "新北市中和區", "臺中市中區", "臺北市萬華區", "臺北市中正區",
        "臺中市西區", "高雄市旗津區", "臺北市信義區", "臺中市南區", "高雄市三民區", 
        "高雄市鹽埕區")

m.data = motor[(motor[, 4] %in% con) == TRUE, ]
View(m.data)
c.data = car[(car[, 4] %in% con) == TRUE, ]
rownames(m.data) = NULL
rownames(c.data) = NULL
unique(c.data[, 4])
write.csv(m.data, "Documents/real data/usage/mot.csv", fileEncoding="big-5")
write.csv(c.data, "Documents/real data/usage/car.csv", fileEncoding="big-5")
###################################################################################################

m = read.csv("Documents/real data/usage/mot.csv", fileEncoding="big-5")
c = read.csv("Documents/real data/usage/car.csv", fileEncoding="big-5")
head(m)
unique(m[, 3])

tab.m = table(m$County, m$Month, m$Year)
tab.c = table(c$County, c$Month, c$Year)

tmp.m = c()
tmp.c = c()
tmp.sq = c()
sq = c(101:110)
for (i in 1:10){
  tmp.sq = c(tmp.sq, rep(sq[i], 21))
  tmp.m = rbind(tmp.m, as.matrix(tab.m[, , i]))
  tmp.c = rbind(tmp.c, as.matrix(tab.c[, , i]))  
}

final1 = cbind(tmp.sq, tmp.m, tmp.c)
colnames(final1) = c("Year", "N11", "N12", "N13", "N14",
                     "N15", "N16", "N17", "N18",
                     "N19", "N110", "N111", "N112",
                     "N21", "N22", "N23", "N24",
                     "N25", "N26", "N27", "N28",
                     "N29", "N210", "N211", "N212")
write.csv(final1, "Documents/real data/Final/final_file.csv", fileEncoding="big-5")
