lowlander_fixture <- function() {
  list(
    latin_mod = matrix(c(0, 10, 1, 11, 2, 12, 3, 13), ncol = 2, byrow = TRUE),
    keep = 1:4,
    output = c(4, 1, 3, 2),
    best = c(first = 1, second = 11),
    lower = c(0, 10),
    upper = c(3, 13),
    liketype = "min"
  )
}
