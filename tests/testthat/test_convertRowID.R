context("convertRowID")

#
# Case1 : 1 vs 1
#
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "B", "C", "D")
LtoR <- rbind(
  c("A", "3"),
  c("B", "2"),
  c("C", "4"),
  c("D", "7"))

out <- convertRowID(input, rowID, LtoR, "sum")
expect_true(identical(rownames(out$output), LtoR[,2]))

out <- convertRowID(input, rowID, LtoR, "mean")
expect_true(identical(rownames(out$output), LtoR[,2]))

out <- convertRowID(input, rowID, LtoR, "large.mean")
expect_true(identical(rownames(out$output), LtoR[,2]))

out <- convertRowID(input, rowID, LtoR, "large.var")
expect_true(identical(rownames(out$output), LtoR[,2]))

out <- convertRowID(input, rowID, LtoR, "large.cv2")
expect_true(identical(rownames(out$output), LtoR[,2]))

# Answer
#   [,1] [,2] [,3] [,4] [,5]
# 3    1    5    9   13   17
# 2    2    6   10   14   18
# 4    3    7   11   15   19
# 7    4    8   12   16   20


# Case2 : 1 vs Many
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "B", "C", "D")
LtoR <- rbind(
  c("A", "3"),
  c("B", "2"), # If cannot decide which element should select, the first one is selected
  c("B", "6"),
  c("C", "4"),
  c("D", "7"))

out <- convertRowID(input, rowID, LtoR, "sum")
expect_true(identical(rownames(out$output), c("3", "2", "4", "7")))

out <- convertRowID(input, rowID, LtoR, "mean")
expect_true(identical(rownames(out$output), c("3", "2", "4", "7")))

out <- convertRowID(input, rowID, LtoR, "large.mean")
expect_true(identical(rownames(out$output), c("3", "2", "4", "7")))

out <- convertRowID(input, rowID, LtoR, "large.var")
expect_true(identical(rownames(out$output), c("3", "2", "4", "7")))

out <- convertRowID(input, rowID, LtoR, "large.cv2")
expect_true(identical(rownames(out$output), c("3", "2", "4", "7")))

# Answer
#   [,1] [,2] [,3] [,4] [,5]
# 3    1    5    9   13   17
# 2    2    6   10   14   18
# 4    3    7   11   15   19
# 7    4    8   12   16   20


# Case3 : Many vs 1
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "A", "B", "D")
LtoR <- rbind(
  c("A", "3"),
  c("B", "2"),
  c("C", "4"),
  c("D", "7"))

out <- convertRowID(input, rowID, LtoR, "sum")
expect_true(identical(rownames(out$output), c("3", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "mean")
expect_true(identical(rownames(out$output), c("3", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "large.mean")
expect_true(identical(rownames(out$output), c("3", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "large.var")
expect_true(identical(rownames(out$output), c("3", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "large.cv2")
expect_true(identical(rownames(out$output), c("3", "2", "7")))

# Answer
#   [,1] [,2] [,3] [,4] [,5]
# 3    1    5    9   13   17
# 2    3    7   11   15   19
# 7    4    8   12   16   20


# Case4 : Many vs Many
input <- matrix(1:20, nrow=4, ncol=5)
rowID <- c("A", "A", "B", "D")
LtoR <- rbind(
  c("A", "5"),
  c("B", "2"), # If cannot decide which element should select, the first one is selected
  c("B", "6"),
  c("C", "4"),
  c("D", "7"))

out <- convertRowID(input, rowID, LtoR, "sum")
expect_true(identical(rownames(out$output), c("5", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "mean")
expect_true(identical(rownames(out$output), c("5", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "large.mean")
expect_true(identical(rownames(out$output), c("5", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "large.var")
expect_true(identical(rownames(out$output), c("5", "2", "7")))

out <- convertRowID(input, rowID, LtoR, "large.cv2")
expect_true(identical(rownames(out$output), c("5", "2", "7")))

# Answer
#   [,1] [,2] [,3] [,4] [,5]
# 5    1    5    9   13   17
# 2    3    7   11   15   19
# 7    4    8   12   16   20
