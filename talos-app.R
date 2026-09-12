 

required_packages <- c(
  "shiny",
  "jsonlite",
  "dplyr",
  "tidyr",
  "ggplot2",
  "plotly",
  "DT",
  "stringr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(shiny)
library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(stringr)

options(shiny.maxRequestSize = 4096 * 1024^2)



# 1.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

safe_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x)[1])
}


safe_num_vec <- function(x) {
  if (is.null(x)) return(numeric(0))
  suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
}


safe_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x)[1]
}

rms <- function(x) {
  sqrt(mean(x^2, na.rm = TRUE))
}

cv_percent <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2) return(NA_real_)

  m <- mean(x)

  if (!is.finite(m) || m == 0) return(NA_real_)

  100 * sd(x) / abs(m)
}

asymmetry_index <- function(left, right) {
  if (!is.finite(left) || !is.finite(right)) return(NA_real_)

  denom <- (abs(left) + abs(right)) / 2

  if (denom == 0) return(NA_real_)

  100 * abs(left - right) / denom
}

fmt <- function(x, digits = 2, suffix = "") {
  if (length(x) == 0 || !is.finite(x)) return("NA")

  paste0(
    format(
      round(x, digits),
      nsmall = digits,
      trim = TRUE
    ),
    suffix
  )
}

clean_choice <- function(x, fallback = "Unknown") {
  x <- as.character(x)

  x[
    is.na(x) |
      trimws(x) == ""
  ] <- fallback

  x
}

pathology_display <- function(group, pathology_key, pathology) {

  result <- pathology_key

  missing_key <- is.na(result) | trimws(result) == ""

  result[missing_key] <- pathology[missing_key]

  missing_result <- is.na(result) | trimws(result) == ""

  result[
    missing_result &
      !is.na(group) &
      tolower(group) == "healthy"
  ] <- "Healthy"

  result[
    is.na(result) |
      trimws(result) == ""
  ] <- "Unknown"

  result
}


# 3.

read_metadata <- function(path) {

  jsonlite::fromJSON(
    path,
    simplifyVector = TRUE
  )
}

metadata_row <- function(path, trial_id) {

  m <- read_metadata(path)

  tibble(
    trial_id = trial_id,
    subject = safe_chr(m$subject),
    age = safe_num(m$age),
    gender = safe_chr(m$gender),
    height = safe_num(m$height),
    weight = safe_num(m$weight),
    BMI = safe_num(m$BMI),
    group = safe_chr(m$group),
    pathology = safe_chr(m$pathology),
    pathology_key = safe_chr(m$pathologyKey),
    laterality = safe_chr(m$laterality),
    clinical_deficit_side = safe_chr(m$clinicalDeficitSide),
    evaluation_name = safe_chr(m$evaluationScoreName),
    evaluation_value = safe_num(m$evaluationScoreValue),
    session = safe_num(m$session),
    trial = safe_num(m$trial),
    frequency = safe_num(m$freq),
    protocol = safe_chr(m$protocol)
  )
}


# 4. 

events_to_df <- function(x, side, fs) {

  if (is.null(x) || length(x) == 0) {
    return(
      tibble(
        side = character(),
        event_number = integer(),
        swing_start = numeric(),
        swing_end = numeric(),
        swing_start_time = numeric(),
        swing_end_time = numeric()
      )
    )
  }

  if (is.data.frame(x) || is.matrix(x)) {
    m <- as.matrix(x)
  } else {
    m <- matrix(
      as.numeric(unlist(x)),
      ncol = 2,
      byrow = TRUE
    )
  }

  if (ncol(m) != 2) {
    stop("Gait-event annotations are not in the expected two-column format.")
  }

  tibble(
    side = side,
    event_number = seq_len(nrow(m)),
    swing_start = as.numeric(m[, 1]),
    swing_end = as.numeric(m[, 2]),
    swing_start_time = as.numeric(m[, 1]) / fs,
    swing_end_time = as.numeric(m[, 2]) / fs
  )
}

classify_event_phase <- function(events, uturn_start, uturn_end) {

  events |>
    mutate(
      phase = case_when(
        swing_start < uturn_start &
          swing_end < uturn_start ~ "Straight — Go",

        swing_start > uturn_end &
          swing_end > uturn_end ~ "Straight — Back",

        TRUE ~ "U-turn"
      )
    )
}

calculate_side_cycles <- function(events, fs) {

  straight_events <- events |>
    filter(
      phase %in% c(
        "Straight — Go",
        "Straight — Back"
      )
    )

  if (nrow(straight_events) == 0) {
    return(tibble())
  }

  pieces <- list()

  for (
    phase_name in c(
      "Straight — Go",
      "Straight — Back"
    )
  ) {

    e <- straight_events |>
      filter(
        phase == phase_name
      ) |>
      arrange(
        swing_start
      )

    if (nrow(e) == 0) next

    pieces[[phase_name]] <- e |>
      mutate(
        cycle = row_number(),
        swing_time = (swing_end - swing_start) / fs,
        next_swing_start = lead(swing_start),
        stance_time = (next_swing_start - swing_end) / fs,
        stride_time = (next_swing_start - swing_start) / fs
      ) |>
      select(
        side,
        phase,
        cycle,
        event_number,
        swing_start,
        swing_end,
        swing_time,
        stance_time,
        stride_time
      )
  }

  bind_rows(pieces)
}

summarize_side_cycles <- function(cycles) {

  if (nrow(cycles) == 0) {
    return(
      tibble(
        n_swings = 0,
        n_complete_strides = 0,
        mean_swing_s = NA_real_,
        mean_stance_s = NA_real_,
        mean_stride_s = NA_real_,
        stride_cv_percent = NA_real_,
        swing_cv_percent = NA_real_,
        stance_cv_percent = NA_real_
      )
    )
  }

  tibble(
    n_swings = sum(
      is.finite(
        cycles$swing_time
      )
    ),
    n_complete_strides = sum(
      is.finite(
        cycles$stride_time
      )
    ),
    mean_swing_s = mean(
      cycles$swing_time,
      na.rm = TRUE
    ),
    mean_stance_s = mean(
      cycles$stance_time,
      na.rm = TRUE
    ),
    mean_stride_s = mean(
      cycles$stride_time,
      na.rm = TRUE
    ),
    stride_cv_percent = cv_percent(
      cycles$stride_time
    ),
    swing_cv_percent = cv_percent(
      cycles$swing_time
    ),
    stance_cv_percent = cv_percent(
      cycles$stance_time
    )
  )
}


# 5. 

natural_file_order <- function(x) {
  key <- stringr::str_replace_all(
    tolower(x),
    "\\d+",
    function(z) sprintf("%010d", as.integer(z))
  )
  order(key)
}

read_processed <- function(path, fs, display_names = NULL) {

  if (length(path) < 1) {
    stop("No processed IMU file was supplied.")
  }

  if (is.null(display_names)) {
    display_names <- basename(path)
  }

  if (length(path) != length(display_names)) {
    stop("Processed file names and paths do not match.")
  }

  ord <- natural_file_order(display_names)
  path <- path[ord]
  display_names <- display_names[ord]

  pieces <- lapply(seq_along(path), function(i) {
    dat <- read.delim(
      path[i],
      header = TRUE,
      sep = "\t",
      check.names = FALSE
    )

    dat$.source_file <- display_names[i]
    dat$.part_number <- i
    dat
  })

  schemas <- lapply(pieces, function(x) setdiff(names(x), c(".source_file", ".part_number")))
  reference_schema <- schemas[[1]]

  bad_schema <- which(!vapply(
    schemas,
    function(x) identical(x, reference_schema),
    logical(1)
  ))

  if (length(bad_schema) > 0) {
    stop(
      paste0(
        "Processed IMU part ",
        bad_schema[1],
        " has different columns from the first part."
      )
    )
  }

  dat <- dplyr::bind_rows(pieces)

  required_columns <- c(
    "PacketCounter",
    "HE_FreeAcc_X","HE_FreeAcc_Y","HE_FreeAcc_Z",
    "HE_Gyr_X","HE_Gyr_Y","HE_Gyr_Z",
    "LB_FreeAcc_X","LB_FreeAcc_Y","LB_FreeAcc_Z",
    "LB_Gyr_X","LB_Gyr_Y","LB_Gyr_Z",
    "LF_FreeAcc_X","LF_FreeAcc_Y","LF_FreeAcc_Z",
    "LF_Gyr_X","LF_Gyr_Y","LF_Gyr_Z",
    "RF_FreeAcc_X","RF_FreeAcc_Y","RF_FreeAcc_Z",
    "RF_Gyr_X","RF_Gyr_Y","RF_Gyr_Z"
  )

  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0) {
    stop(
      paste(
        "Processed file is missing:",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  original_counter <- suppressWarnings(as.numeric(dat$PacketCounter))

  counter_ok <- (
    all(is.finite(original_counter)) &&
    
    length(original_counter) > 1 &&
    all(diff(original_counter) > 0)
  )

  if (!counter_ok) {
    dat$PacketCounter <- seq.int(0, nrow(dat) - 1)
  }

  dat |>
    mutate(
      sample = as.numeric(PacketCounter),
      time = sample / fs,

      HE_FreeAcc_mag = sqrt(HE_FreeAcc_X^2 + HE_FreeAcc_Y^2 + HE_FreeAcc_Z^2),
      LB_FreeAcc_mag = sqrt(LB_FreeAcc_X^2 + LB_FreeAcc_Y^2 + LB_FreeAcc_Z^2),
      LF_FreeAcc_mag = sqrt(LF_FreeAcc_X^2 + LF_FreeAcc_Y^2 + LF_FreeAcc_Z^2),
      RF_FreeAcc_mag = sqrt(RF_FreeAcc_X^2 + RF_FreeAcc_Y^2 + RF_FreeAcc_Z^2),

      HE_Gyr_mag = sqrt(HE_Gyr_X^2 + HE_Gyr_Y^2 + HE_Gyr_Z^2),
      LB_Gyr_mag = sqrt(LB_Gyr_X^2 + LB_Gyr_Y^2 + LB_Gyr_Z^2),
      LF_Gyr_mag = sqrt(LF_Gyr_X^2 + LF_Gyr_Y^2 + LF_Gyr_Z^2),
      RF_Gyr_mag = sqrt(RF_Gyr_X^2 + RF_Gyr_Y^2 + RF_Gyr_Z^2)
    )
}

get_walking_windows <- function(
  left_events,
  right_events,
  uturn_start,
  uturn_end
) {

  all_events <- bind_rows(
    left_events,
    right_events
  )

  go_events <- all_events |>
    filter(
      phase == "Straight — Go"
    )

  back_events <- all_events |>
    filter(
      phase == "Straight — Back"
    )

  windows <- tibble(
    phase = character(),
    start_sample = numeric(),
    end_sample = numeric()
  )

  if (nrow(go_events) > 0) {

    windows <- bind_rows(
      windows,
      tibble(
        phase = "Straight — Go",
        start_sample = min(
          go_events$swing_start,
          na.rm = TRUE
        ),
        end_sample = uturn_start - 1
      )
    )
  }

  if (nrow(back_events) > 0) {

    windows <- bind_rows(
      windows,
      tibble(
        phase = "Straight — Back",
        start_sample = uturn_end + 1,
        end_sample = max(
          back_events$swing_end,
          na.rm = TRUE
        )
      )
    )
  }

  windows
}

filter_straight_samples <- function(dat, windows) {

  if (nrow(windows) == 0) {
    return(
      dat[0, , drop = FALSE]
    )
  }

  keep <- rep(
    FALSE,
    nrow(dat)
  )

  for (i in seq_len(nrow(windows))) {

    keep <- keep |
      (
        dat$sample >= windows$start_sample[i] &
          dat$sample <= windows$end_sample[i]
      )
  }

  dat[
    keep,
    ,
    drop = FALSE
  ]
}

sensor_summary <- function(walking_data) {

  if (nrow(walking_data) == 0) {
    return(tibble())
  }

  left_acc <- rms(
    walking_data$LF_FreeAcc_mag
  )

  right_acc <- rms(
    walking_data$RF_FreeAcc_mag
  )

  left_gyr <- rms(
    walking_data$LF_Gyr_mag
  )

  right_gyr <- rms(
    walking_data$RF_Gyr_mag
  )

  trunk_acc <- rms(
    walking_data$LB_FreeAcc_mag
  )

  head_acc <- rms(
    walking_data$HE_FreeAcc_mag
  )

  foot_mean <- mean(
    c(
      left_acc,
      right_acc
    ),
    na.rm = TRUE
  )

  left_right_coupling <- abs(
    cor(
      walking_data$LF_FreeAcc_mag,
      walking_data$RF_FreeAcc_mag,
      use = "complete.obs"
    )
  )

  foot_trunk_coupling <- mean(
    c(
      abs(
        cor(
          walking_data$LF_FreeAcc_mag,
          walking_data$LB_FreeAcc_mag,
          use = "complete.obs"
        )
      ),
      abs(
        cor(
          walking_data$RF_FreeAcc_mag,
          walking_data$LB_FreeAcc_mag,
          use = "complete.obs"
        )
      )
    ),
    na.rm = TRUE
  )

  tibble(
    Metric = c(
      "Left foot free-acc RMS",
      "Right foot free-acc RMS",
      "Foot free-acc asymmetry",
      "Left foot gyro RMS",
      "Right foot gyro RMS",
      "Foot gyro asymmetry",
      "Trunk free-acc RMS",
      "Head free-acc RMS",
      "Trunk / foot acceleration ratio",
      "Head / trunk acceleration ratio",
      "Left-right foot movement coupling",
      "Mean foot-trunk movement coupling"
    ),

    Value = c(
      left_acc,
      right_acc,
      asymmetry_index(
        left_acc,
        right_acc
      ),
      left_gyr,
      right_gyr,
      asymmetry_index(
        left_gyr,
        right_gyr
      ),
      trunk_acc,
      head_acc,
      trunk_acc / foot_mean,
      head_acc / trunk_acc,
      left_right_coupling,
      foot_trunk_coupling
    )
  )
}


# 6. 

analyze_trial <- function(
  processed_path,
  metadata_path,
  processed_names = NULL
) {

  m <- read_metadata(
    metadata_path
  )

  fs <- safe_num(
    m$freq
  )

  if (
    !is.finite(fs) ||
      fs <= 0
  ) {
    stop(
      "Metadata does not contain a valid sampling frequency."
    )
  }

  dat <- read_processed(
    processed_path,
    fs,
    display_names = processed_names
  )

  uturn <- as.numeric(
    m$uturnBoundaries
  )

  if (length(uturn) < 2) {
    stop(
      "Metadata does not contain valid U-turn boundaries."
    )
  }

  uturn_start <- uturn[1]
  uturn_end <- uturn[2]

  left_events <- events_to_df(
    m$leftGaitEvents,
    "Left",
    fs
  ) |>
    classify_event_phase(
      uturn_start,
      uturn_end
    )

  right_events <- events_to_df(
    m$rightGaitEvents,
    "Right",
    fs
  ) |>
    classify_event_phase(
      uturn_start,
      uturn_end
    )

  left_cycles <- calculate_side_cycles(
    left_events,
    fs
  )

  right_cycles <- calculate_side_cycles(
    right_events,
    fs
  )

  left_summary <- summarize_side_cycles(
    left_cycles
  )

  right_summary <- summarize_side_cycles(
    right_cycles
  )

  timing_summary <- tibble(
    Metric = c(
      "Stride time",
      "Swing time",
      "Stance time",
      "Stride variability"
    ),

    Left = c(
      left_summary$mean_stride_s,
      left_summary$mean_swing_s,
      left_summary$mean_stance_s,
      left_summary$stride_cv_percent
    ),

    Right = c(
      right_summary$mean_stride_s,
      right_summary$mean_swing_s,
      right_summary$mean_stance_s,
      right_summary$stride_cv_percent
    ),

    Asymmetry = c(
      asymmetry_index(
        left_summary$mean_stride_s,
        right_summary$mean_stride_s
      ),

      asymmetry_index(
        left_summary$mean_swing_s,
        right_summary$mean_swing_s
      ),

      asymmetry_index(
        left_summary$mean_stance_s,
        right_summary$mean_stance_s
      ),

      asymmetry_index(
        left_summary$stride_cv_percent,
        right_summary$stride_cv_percent
      )
    )
  )

  windows <- get_walking_windows(
    left_events,
    right_events,
    uturn_start,
    uturn_end
  )

  walking_data <- filter_straight_samples(
    dat,
    windows
  )

  sensors <- sensor_summary(
    walking_data
  )

  cycles <- bind_rows(
    left_cycles,
    right_cycles
  ) |>
    arrange(
      phase,
      side,
      cycle
    )

  list(
    metadata = m,
    data = dat,
    walking_data = walking_data,
    left_events = left_events,
    right_events = right_events,
    left_cycles = left_cycles,
    right_cycles = right_cycles,
    cycles = cycles,
    left_summary = left_summary,
    right_summary = right_summary,
    timing_summary = timing_summary,
    sensor_summary = sensors,
    windows = windows,
    fs = fs,
    uturn_start = uturn_start,
    uturn_end = uturn_end
  )
}



# 7. 


best_gyro_axis <- function(dat, prefix) {

  candidates <- paste0(
    prefix,
    "_Gyr_",
    c(
      "X",
      "Y",
      "Z"
    )
  )

  sds <- vapply(
    candidates,
    function(x) {
      sd(
        dat[[x]],
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  candidates[
    which.max(sds)
  ]
}

build_intervals <- function(events) {

  if (nrow(events) == 0) {
    return(tibble())
  }

  pieces <- list()

  for (i in seq_len(nrow(events))) {

    pieces[[
      length(pieces) + 1
    ]] <- tibble(
      state = "Swing",
      xmin = events$swing_start_time[i],
      xmax = events$swing_end_time[i]
    )

    if (i < nrow(events)) {

      if (
        events$phase[i] ==
          events$phase[i + 1]
      ) {

        pieces[[
          length(pieces) + 1
        ]] <- tibble(
          state = "Stance",
          xmin = events$swing_end_time[i],
          xmax = events$swing_start_time[i + 1]
        )
      }
    }
  }

  bind_rows(
    pieces
  )
}

make_gait_plot <- function(result) {

  dat <- result$data

  left_signal <- best_gyro_axis(dat, "LF")
  right_signal <- best_gyro_axis(dat, "RF")

  left_df <- data.frame(
    time = dat$time,
    angular_velocity = dat[[left_signal]],
    side = "Left foot",
    stringsAsFactors = FALSE
  )

  right_df <- data.frame(
    time = dat$time,
    angular_velocity = dat[[right_signal]],
    side = "Right foot",
    stringsAsFactors = FALSE
  )

  plot_df <- dplyr::bind_rows(left_df, right_df)
  plot_df$side <- factor(
    plot_df$side,
    levels = c("Left foot", "Right foot")
  )

  make_side_intervals <- function(events, side_label) {
    x <- build_intervals(events)
    if (nrow(x) == 0) return(x)
    x$side <- side_label
    x
  }

  intervals <- dplyr::bind_rows(
    make_side_intervals(result$left_events, "Left foot"),
    make_side_intervals(result$right_events, "Right foot")
  )

  event_lines <- dplyr::bind_rows(
    data.frame(
      x = c(
        result$left_events$swing_start_time,
        result$left_events$swing_end_time
      ),
      side = "Left foot"
    ),
    data.frame(
      x = c(
        result$right_events$swing_start_time,
        result$right_events$swing_end_time
      ),
      side = "Right foot"
    )
  )

  event_lines$side <- factor(
    event_lines$side,
    levels = c("Left foot", "Right foot")
  )

  turn_df <- data.frame(
    xmin = result$uturn_start / result$fs,
    xmax = result$uturn_end / result$fs,
    side = factor(
      c("Left foot", "Right foot"),
      levels = c("Left foot", "Right foot")
    )
  )

  p <- ggplot(
    plot_df,
    aes(
      x = time,
      y = angular_velocity
    )
  ) +

    geom_rect(
      data = turn_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "#a8793d",
      alpha = 0.10
    ) +

    geom_rect(
      data = intervals,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = state
      ),
      inherit.aes = FALSE,
      alpha = 0.12
    ) +

    geom_line(
      colour = "#173047",
      linewidth = 0.45
    ) +

    geom_vline(
      data = event_lines,
      aes(xintercept = x),
      inherit.aes = FALSE,
      linetype = "dashed",
      colour = "#6f6255",
      linewidth = 0.28,
      alpha = 0.70
    ) +

    facet_wrap(
      ~side,
      ncol = 1,
      scales = "free_y"
    ) +

    scale_fill_manual(
      values = c(
        "Swing" = "#d9a0a0",
        "Stance" = "#9ab99b"
      ),
      guide = "none"
    ) +

    labs(
      x = "Time (s)",
      y = "Angular velocity (rad/s)",
      fill = NULL,
      caption = paste0(
        "Left foot: ", left_signal,
        "  •  Right foot: ", right_signal,
        "  •  Bronze shading: U-turn"
      )
    ) +

    theme_minimal(
      base_size = 12
    ) +

    theme(
      strip.text = element_text(
        family = "serif",
        face = "bold",
        size = 12,
        colour = "#173047"
      ),
      strip.background = element_rect(
        fill = "#f3eadc",
        colour = "#d9c7ab",
        linewidth = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "#e8e0d4",
        linewidth = 0.35
      ),
      legend.position = "none",
      plot.caption = element_text(
        hjust = 0,
        colour = "#80694d",
        size = 9
      )
    )

  ggplotly(
    p,
    tooltip = c("x", "y")
  )
}

make_sensor_plot <- function(result) {

  dat <- result$data
  fs <- result$fs
  meta <- result$metadata

  sensor_map <- c(
    HE = "Head",
    LB = "Lower back",
    LF = "Left foot",
    RF = "Right foot"
  )

  
  find_xyz <- function(prefix, type = c("FreeAcc", "Acc")) {
    type <- match.arg(type)
    candidates <- list(
      paste0(prefix, "_", type, "_", c("X", "Y", "Z")),
      paste0(prefix, ".", type, ".", c("X", "Y", "Z")),
      paste0(prefix, "_", type, ".", c("X", "Y", "Z")),
      paste0(prefix, ".", type, "_", c("X", "Y", "Z"))
    )

    for (cols in candidates) {
      if (all(cols %in% names(dat))) return(cols)
    }

    NULL
  }

  plot_rows <- lapply(names(sensor_map), function(prefix) {

    cols <- find_xyz(prefix, "FreeAcc")
    signal_type <- "Free acceleration"

    if (is.null(cols)) {
      cols <- find_xyz(prefix, "Acc")
      signal_type <- "Acceleration"
    }

    if (is.null(cols)) return(NULL)

    x <- suppressWarnings(as.numeric(dat[[cols[1]]]))
    y <- suppressWarnings(as.numeric(dat[[cols[2]]]))
    z <- suppressWarnings(as.numeric(dat[[cols[3]]]))

    mag <- sqrt(x^2 + y^2 + z^2)

    data.frame(
      time_s = seq_along(mag) / fs,
      magnitude = mag,
      sensor = sensor_map[[prefix]],
      signal_type = signal_type,
      stringsAsFactors = FALSE
    )
  })

  plot_df <- dplyr::bind_rows(plot_rows)

  if (nrow(plot_df) == 0) {
    stop(
      paste0(
        "No compatible acceleration columns were found in the analyzed data. ",
        "Expected names include LF_FreeAcc_X/Y/Z or LF_Acc_X/Y/Z. ",
        "Available columns begin: ",
        paste(utils::head(names(dat), 12), collapse = ", ")
      )
    )
  }

  plot_df$sensor <- factor(
    plot_df$sensor,
    levels = c("Head", "Lower back", "Left foot", "Right foot")
  )

  turn <- safe_num_vec(meta$uturnBoundaries)
  turn_df <- NULL

  if (length(turn) >= 2 && all(is.finite(turn[1:2]))) {
    turn_df <- data.frame(
      xmin = turn[1] / fs,
      xmax = turn[2] / fs,
      ymin = -Inf,
      ymax = Inf
    )
  }

  p <- ggplot(plot_df, aes(x = time_s, y = magnitude))

  if (!is.null(turn_df)) {
    p <- p + geom_rect(
      data = turn_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#a8793d",
      alpha = 0.10
    )
  }

  signal_label <- if (all(plot_df$signal_type == "Free acceleration")) {
    "Free-acceleration magnitude"
  } else if (all(plot_df$signal_type == "Acceleration")) {
    "Acceleration magnitude"
  } else {
    "Acceleration magnitude"
  }

  p +
    geom_line(linewidth = 0.45, colour = "#173047") +
    facet_wrap(~sensor, ncol = 1, scales = "free_y") +
    labs(
      x = "Time (s)",
      y = signal_label,
      caption = "Bronze shading marks the annotated U-turn."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(
        family = "serif",
        face = "bold",
        size = 12,
        colour = "#173047"
      ),
      strip.background = element_rect(
        fill = "#f3eadc",
        colour = "#d9c7ab",
        linewidth = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "#e8e0d4",
        linewidth = 0.35
      ),
      axis.title = element_text(colour = "#31495c"),
      axis.text = element_text(colour = "#566574"),
      plot.caption = element_text(
        hjust = 0,
        colour = "#80694d",
        size = 9
      )
    )
}




# 8. 

fmt_value <- function(x, digits = 2, suffix = "") {
  x <- suppressWarnings(as.numeric(x)[1])
  if (!is.finite(x)) return("—")
  paste0(format(round(x, digits), nsmall = digits, trim = TRUE), suffix)
}

get_metric <- function(tbl, metric, column) {
  z <- tbl[tbl$Metric == metric, column, drop = TRUE]
  if (length(z) == 0) return(NA_real_)
  suppressWarnings(as.numeric(z[1]))
}


make_interpretation <- function(r) {
  t <- r$timing_summary

  stride_asym <- get_metric(t, "Stride time", "Asymmetry")
  swing_asym  <- get_metric(t, "Swing time", "Asymmetry")
  stance_asym <- get_metric(t, "Stance time", "Asymmetry")

  left_swing   <- get_metric(t, "Swing time", "Left")
  right_swing  <- get_metric(t, "Swing time", "Right")
  left_stance  <- get_metric(t, "Stance time", "Left")
  right_stance <- get_metric(t, "Stance time", "Right")

  notes <- character(0)

  if (all(is.finite(c(stride_asym, swing_asym, stance_asym)))) {
    notes <- c(
      notes,
      sprintf(
        "Stride timing is nearly matched between sides (%.1f%% difference); the larger differences are within the stride cycle.",
        stride_asym
      ),
      sprintf(
        "Swing: %.1f%% difference • Stance: %.1f%% difference.",
        swing_asym, stance_asym
      )
    )
  }

  if (all(is.finite(c(left_swing, right_swing)))) {
    notes <- c(
      notes,
      sprintf(
        "Swing is longer on the %s (L %.3f s • R %.3f s).",
        ifelse(left_swing > right_swing, "right", "left"),
        left_swing, right_swing
      )
    )
  }

  if (all(is.finite(c(left_stance, right_stance)))) {
    notes <- c(
      notes,
      sprintf(
        "Stance is longer on the %s (L %.3f s • R %.3f s).",
        ifelse(left_stance > right_stance, "left", "right"),
        left_stance, right_stance
      )
    )
  }

  notes
}

combined_summary_table <- function(r) {
  timing <- r$timing_summary
  timing_out <- data.frame(
    Section = "Gait timing",
    Metric = timing$Metric,
    Left = timing$Left,
    Right = timing$Right,
    Value = NA_real_,
    `Asymmetry (%)` = timing$Asymmetry,
    check.names = FALSE
  )

  sensor <- r$sensor_summary
  sensor_out <- data.frame(
    Section = "Movement",
    Metric = sensor$Metric,
    Left = NA_real_,
    Right = NA_real_,
    Value = sensor$Value,
    `Asymmetry (%)` = NA_real_,
    check.names = FALSE
  )

  rbind(timing_out, sensor_out)
}

metric_card <- function(label, value, note = "") {
  div(
    class = "metric-card",
    div(class = "metric-label", label),
    div(class = "metric-value", value),
    div(class = "metric-note", note)
  )
}




TALOS_VERSION <- "1.0.0-stable-no-variability-plot"

TALOS_REQUIRED_META <- c(
  "subject", "group", "pathology", "freq",
  "uturnBoundaries", "leftGaitEvents", "rightGaitEvents"
)

TALOS_SENSOR_PREFIXES <- c("HE", "LB", "LF", "RF")

TALOS_SENSOR_LABELS <- c(
  HE = "Head",
  LB = "Lower back",
  LF = "Left foot",
  RF = "Right foot"
)

talos_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

talos_trim_names <- function(x) {
  x <- trimws(x)
  x <- gsub("\ufeff", "", x, fixed = TRUE)
  x
}

talos_standardize_names <- function(nms) {
  nms <- talos_trim_names(nms)
  nms <- gsub("\\.", "_", nms)
  nms <- gsub("__+", "_", nms)
  nms
}

talos_read_text_head <- function(path, n = 8) {
  tryCatch(
    readLines(path, n = n, warn = FALSE, encoding = "UTF-8"),
    error = function(e) character(0)
  )
}

talos_detect_separator <- function(path) {
  lines <- talos_read_text_head(path, n = 5)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) return("\t")

  candidate <- lines[1]
  n_tab <- lengths(regmatches(candidate, gregexpr("\t", candidate, fixed = TRUE)))
  n_comma <- lengths(regmatches(candidate, gregexpr(",", candidate, fixed = TRUE)))
  n_semicolon <- lengths(regmatches(candidate, gregexpr(";", candidate, fixed = TRUE)))

  if (n_tab >= n_comma && n_tab >= n_semicolon) return("\t")
  if (n_comma >= n_semicolon) return(",")
  ";"
}

talos_read_table_robust <- function(path) {
  sep <- talos_detect_separator(path)

  dat <- tryCatch(
    read.table(
      path,
      header = TRUE,
      sep = sep,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    ),
    error = function(e) {
      stop(
        paste0("Could not read processed IMU file '", basename(path), "': ", conditionMessage(e)),
        call. = FALSE
      )
    }
  )

  names(dat) <- talos_standardize_names(names(dat))
  dat
}

talos_normalize_metadata <- function(meta) {
  aliases <- list(
    subject = c("subject", "participant", "participant_id", "subject_id"),
    freq = c("freq", "frequency", "sampling_rate", "samplingRate"),
    leftGaitEvents = c("leftGaitEvents", "left_gait_events", "leftEvents"),
    rightGaitEvents = c("rightGaitEvents", "right_gait_events", "rightEvents"),
    uturnBoundaries = c("uturnBoundaries", "uturn_boundaries", "uTurnBoundaries"),
    clinicalDeficitSide = c("clinicalDeficitSide", "clinical_deficit_side", "deficitSide")
  )

  for (target in names(aliases)) {
    if (!is.null(meta[[target]])) next
    for (candidate in aliases[[target]]) {
      if (!is.null(meta[[candidate]])) {
        meta[[target]] <- meta[[candidate]]
        break
      }
    }
  }

  meta
}

talos_validate_metadata <- function(meta) {
  warnings <- character(0)
  errors <- character(0)

  if (!is.list(meta)) {
    return(list(ok = FALSE, errors = "Metadata is not a JSON object.", warnings = warnings))
  }

  meta <- talos_normalize_metadata(meta)

  missing <- TALOS_REQUIRED_META[vapply(
    TALOS_REQUIRED_META,
    function(x) is.null(meta[[x]]),
    logical(1)
  )]

  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing required metadata fields:", paste(missing, collapse = ", ")))
  }

  fs <- suppressWarnings(as.numeric(meta$freq))
  if (length(fs) != 1 || !is.finite(fs) || fs <= 0) {
    errors <- c(errors, "Sampling frequency must be one positive numeric value.")
  } else if (fs < 20 || fs > 2000) {
    warnings <- c(warnings, paste0("Sampling frequency ", fs, " Hz is unusual; verify metadata."))
  }

  turn <- suppressWarnings(as.numeric(unlist(meta$uturnBoundaries)))
  if (length(turn) < 2 || any(!is.finite(turn[1:2]))) {
    errors <- c(errors, "U-turn boundaries must contain two finite sample indices.")
  } else if (turn[2] <= turn[1]) {
    errors <- c(errors, "U-turn end must occur after U-turn start.")
  }

  validate_events <- function(x, side) {
    err <- character(0)
    warn <- character(0)

    if (is.null(x)) {
      return(list(errors = paste(side, "gait events are missing."), warnings = warn))
    }

    mat <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (is.null(mat) || ncol(mat) != 2) {
      return(list(
        errors = paste(side, "gait events must be two-column swing start/end pairs."),
        warnings = warn
      ))
    }

    vals <- suppressWarnings(matrix(as.numeric(mat), ncol = 2))

    if (any(!is.finite(vals))) {
      err <- c(err, paste(side, "gait events contain non-numeric values."))
    }

    if (any(vals[, 2] <= vals[, 1], na.rm = TRUE)) {
      err <- c(err, paste(side, "gait events contain end <= start."))
    }

    if (nrow(vals) > 1 && any(diff(vals[, 1]) <= 0, na.rm = TRUE)) {
      warn <- c(warn, paste(side, "gait events are not strictly increasing."))
    }

    list(errors = err, warnings = warn)
  }

  l <- validate_events(meta$leftGaitEvents, "Left")
  r <- validate_events(meta$rightGaitEvents, "Right")

  errors <- c(errors, l$errors, r$errors)
  warnings <- c(warnings, l$warnings, r$warnings)

  list(
    ok = length(errors) == 0,
    errors = unique(errors),
    warnings = unique(warnings),
    metadata = meta
  )
}

talos_find_sensor_xyz <- function(dat, prefix, signal_type = c("FreeAcc", "Acc", "Gyr")) {
  signal_type <- match.arg(signal_type)
  wanted <- paste0(prefix, "_", signal_type, "_", c("X", "Y", "Z"))
  current <- names(dat)
  idx <- match(tolower(wanted), tolower(current))
  if (all(!is.na(idx))) current[idx] else NULL
}

talos_signal_inventory <- function(dat) {
  dplyr::bind_rows(lapply(TALOS_SENSOR_PREFIXES, function(prefix) {
    data.frame(
      sensor = TALOS_SENSOR_LABELS[[prefix]],
      prefix = prefix,
      free_acceleration = !is.null(talos_find_sensor_xyz(dat, prefix, "FreeAcc")),
      acceleration = !is.null(talos_find_sensor_xyz(dat, prefix, "Acc")),
      gyroscope = !is.null(talos_find_sensor_xyz(dat, prefix, "Gyr")),
      stringsAsFactors = FALSE
    )
  }))
}

talos_coerce_numeric_signals <- function(dat) {
  candidates <- grep(
    "(PacketCounter|_(FreeAcc|Acc|Gyr)_[XYZ])$",
    names(dat),
    value = TRUE
  )

  for (nm in candidates) {
    if (!is.numeric(dat[[nm]])) {
      dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
    }
  }

  dat
}

talos_validate_processed_data <- function(dat, meta = NULL) {
  warnings <- character(0)
  errors <- character(0)

  if (!is.data.frame(dat) || nrow(dat) == 0) {
    return(list(ok = FALSE, errors = "Processed IMU data contains no rows.", warnings = warnings))
  }

  names(dat) <- talos_standardize_names(names(dat))
  dat <- talos_coerce_numeric_signals(dat)

  dups <- unique(names(dat)[duplicated(names(dat))])
  if (length(dups) > 0) {
    errors <- c(errors, paste("Duplicate columns:", paste(dups, collapse = ", ")))
  }

  inventory <- talos_signal_inventory(dat)

  if (!inventory$gyroscope[inventory$prefix == "LF"]) {
    errors <- c(errors, "Left-foot gyroscope XYZ channels are missing.")
  }

  if (!inventory$gyroscope[inventory$prefix == "RF"]) {
    errors <- c(errors, "Right-foot gyroscope XYZ channels are missing.")
  }

  if (!any(inventory$free_acceleration | inventory$acceleration)) {
    warnings <- c(
      warnings,
      "No complete XYZ acceleration signal was found; the segmental acceleration plot may be unavailable."
    )
  }

  numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  if (length(numeric_cols) > 0) {
    bad_rates <- vapply(
      dat[numeric_cols],
      function(x) mean(!is.finite(x)),
      numeric(1)
    )

    bad <- names(bad_rates)[bad_rates > 0.05]

    if (length(bad) > 0) {
      warnings <- c(
        warnings,
        paste0(
          "More than 5% missing/non-finite values in: ",
          paste(utils::head(bad, 8), collapse = ", "),
          if (length(bad) > 8) " …" else ""
        )
      )
    }
  }

  if (!is.null(meta)) {
    turn <- suppressWarnings(as.numeric(unlist(meta$uturnBoundaries)))

    if (length(turn) >= 2 && all(is.finite(turn[1:2]))) {
      if (turn[1] < 0 || turn[2] > nrow(dat) + 5) {
        errors <- c(errors, "U-turn boundaries fall outside the processed data range.")
      }
    }

    ev <- suppressWarnings(as.numeric(unlist(c(meta$leftGaitEvents, meta$rightGaitEvents))))
    if (length(ev) > 0 && any(is.finite(ev))) {
      mx <- max(ev, na.rm = TRUE)
      if (is.finite(mx) && mx > nrow(dat) + 5) {
        errors <- c(errors, "One or more gait-event indices exceed the processed data length.")
      }
    }
  }

  list(
    ok = length(errors) == 0,
    errors = unique(errors),
    warnings = unique(warnings),
    inventory = inventory,
    data = dat
  )
}

talos_validate_multipart <- function(paths, display_names = basename(paths)) {
  warnings <- character(0)
  errors <- character(0)

  if (length(paths) == 0) {
    return(list(ok = FALSE, errors = "No processed files supplied.", warnings = warnings))
  }

  ord <- natural_file_order(display_names)
  paths <- paths[ord]
  display_names <- display_names[ord]

  pieces <- lapply(paths, talos_read_table_robust)

  schemas <- lapply(pieces, names)
  reference <- schemas[[1]]

  mismatch <- which(!vapply(schemas, function(x) identical(x, reference), logical(1)))

  if (length(mismatch) > 0) {
    errors <- c(
      errors,
      paste0(
        "Processed part '", display_names[mismatch[1]],
        "' has a different schema from the first part."
      )
    )
  }

  if ("PacketCounter" %in% reference && length(pieces) > 1) {
    bounds <- lapply(pieces, function(x) {
      pc <- suppressWarnings(as.numeric(x$PacketCounter))
      c(first = pc[1], last = pc[length(pc)])
    })

    for (i in 2:length(bounds)) {
      prev <- bounds[[i - 1]]
      cur <- bounds[[i]]

      if (all(is.finite(c(prev, cur)))) {
        if (cur["first"] <= prev["last"]) {
          warnings <- c(
            warnings,
            paste0(
              "Part ", i,
              " restarts/overlaps PacketCounter; TALOS will rebuild a continuous counter."
            )
          )
        } else if (cur["first"] > prev["last"] + 1) {
          warnings <- c(
            warnings,
            paste0("PacketCounter gap detected between parts ", i - 1, " and ", i, ".")
          )
        }
      }
    }
  }

  list(
    ok = length(errors) == 0,
    errors = unique(errors),
    warnings = unique(warnings),
    order = ord,
    paths = paths,
    names = display_names
  )
}

talos_file_manifest <- function(paths, display_names = basename(paths)) {
  info <- file.info(paths)

  data.frame(
    file = display_names,
    bytes = as.numeric(info$size),
    modified = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
}

talos_build_provenance <- function(meta, processed_paths, processed_names, data_rows, warnings) {
  list(
    talos_version = TALOS_VERSION,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    platform = R.version$platform,
    subject = safe_chr(meta$subject),
    sampling_rate_hz = safe_num(meta$freq),
    processed_file_count = length(processed_paths),
    processed_files = talos_file_manifest(processed_paths, processed_names),
    combined_rows = data_rows,
    warnings = warnings
  )
}

talos_quality_summary <- function(r) {
  inv <- talos_signal_inventory(r$data)

  data.frame(
    check = c(
      "Processed samples",
      "Left gait events",
      "Right gait events",
      "Gyroscope sensor sets",
      "Acceleration sensor sets",
      "Processed files combined"
    ),
    value = c(
      nrow(r$data),
      nrow(r$left_events),
      nrow(r$right_events),
      sum(inv$gyroscope),
      sum(inv$free_acceleration | inv$acceleration),
      r$processed_file_count
    ),
    stringsAsFactors = FALSE
  )
}

talos_internal_self_test <- function() {
  talos_assert(
    abs(asymmetry_index(1, 1)) < 1e-12,
    "Internal test failed: asymmetry identity."
  )

  talos_assert(
    abs(asymmetry_index(1, 2) - 66.6666667) < 1e-5,
    "Internal test failed: asymmetry formula."
  )

  files <- c("part10.txt", "part2.txt", "part1.txt")
  ordered <- files[natural_file_order(files)]

  talos_assert(
    identical(ordered, c("part1.txt", "part2.txt", "part10.txt")),
    "Internal test failed: natural multipart ordering."
  )

  TRUE
}

tryCatch(
  talos_internal_self_test(),
  error = function(e) {
    stop(
      paste0("TALOS startup self-test failed: ", conditionMessage(e)),
      call. = FALSE
    )
  }
)

analyze_trial_hardened <- function(processed_path, metadata_path, processed_names = NULL) {
  if (is.null(processed_names)) processed_names <- basename(processed_path)

  meta <- read_metadata(metadata_path)
  meta_check <- talos_validate_metadata(meta)

  if (!meta_check$ok) {
    stop(paste(meta_check$errors, collapse = " | "), call. = FALSE)
  }

  multi_check <- talos_validate_multipart(processed_path, processed_names)

  if (!multi_check$ok) {
    stop(paste(multi_check$errors, collapse = " | "), call. = FALSE)
  }

  r <- analyze_trial(
    processed_path,
    metadata_path,
    processed_names = processed_names
  )

  data_check <- talos_validate_processed_data(
    r$data,
    meta = meta_check$metadata
  )

  if (!data_check$ok) {
    stop(paste(data_check$errors, collapse = " | "), call. = FALSE)
  }

  r$data <- data_check$data
  r$processed_file_count <- length(processed_path)
  r$processed_file_names <- processed_names
  r$qc_warnings <- unique(c(
    meta_check$warnings,
    multi_check$warnings,
    data_check$warnings
  ))
  r$signal_inventory <- data_check$inventory
  r$quality_summary <- talos_quality_summary(r)
  r$provenance <- talos_build_provenance(
    meta = meta_check$metadata,
    processed_paths = processed_path,
    processed_names = processed_names,
    data_rows = nrow(r$data),
    warnings = r$qc_warnings
  )

  class(r) <- c("talos_analysis", class(r))
  r
}

app_css <- "
:root{--ink:#18232d;--muted:#66727c;--line:#d8dde1;--soft:#f5f6f6;--paper:#fbfbfa;--navy:#163247;--bronze:#9a7444;}
html,body{background:var(--paper);color:var(--ink);}
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;}
.container-fluid{padding:0!important;max-width:none!important;}
.research-header{height:74px;padding:0 34px;display:flex;align-items:center;justify-content:space-between;background:#fff;border-bottom:1px solid var(--line);}
.wordmark{font-family:Georgia,'Times New Roman',serif;font-size:28px;letter-spacing:.12em;color:var(--navy);}
.wordmark-sub{margin-top:2px;font-size:10px;letter-spacing:.13em;text-transform:uppercase;color:var(--muted);}
.header-rule{width:88px;height:3px;background:var(--bronze);margin-top:7px;}
.research-shell{display:grid;grid-template-columns:255px minmax(0,1fr);min-height:calc(100vh - 74px);}
.research-sidebar{background:#f2f4f4;border-right:1px solid var(--line);padding:25px 20px;}
.side-title{font-size:10px;font-weight:700;letter-spacing:.15em;text-transform:uppercase;color:var(--muted);margin-bottom:14px;}
.upload-box{margin-bottom:17px}.control-label{font-size:11px;font-weight:650;color:var(--ink);margin-bottom:6px}
.btn-default,.btn-file{border:1px solid #aeb7bd!important;background:#fff!important;color:var(--navy)!important;border-radius:3px!important;box-shadow:none!important;font-size:11px!important;}
.form-control{border-radius:3px!important;box-shadow:none!important;border-color:#c8ced2!important;font-size:11px!important;}
.status{margin:18px 0;padding:9px 10px;border-left:3px solid var(--line);background:#fff;font-size:11px;line-height:1.45}
.status.ready{border-left-color:#627e69}.status.error{border-left-color:#a65d55}.status.wait{border-left-color:#9b8a72}
.btn-download{width:100%;margin-top:8px!important;padding:9px 10px!important;border:1px solid var(--navy)!important;background:var(--navy)!important;color:#fff!important;border-radius:3px!important;font-size:11px!important;text-align:left!important;}
.small-note{margin-top:20px;padding-top:15px;border-top:1px solid var(--line);font-size:10px;line-height:1.5;color:var(--muted)}
.workspace{padding:28px 38px 54px;max-width:1380px;margin:0 auto}.empty{margin:90px auto;max-width:600px;text-align:center;color:var(--muted)}
.empty h3{font-family:Georgia,'Times New Roman',serif;color:var(--navy);font-weight:normal}
.trial-strip{display:grid;grid-template-columns:1fr .9fr 1.45fr 1fr .8fr .8fr;background:#fff;border:1px solid var(--line);margin-bottom:26px}
.trial-cell{padding:12px 15px;border-right:1px solid var(--line);min-width:0}.trial-cell:last-child{border-right:0}
.trial-label{font-size:9px;letter-spacing:.13em;text-transform:uppercase;color:var(--muted);margin-bottom:4px}
.trial-value{font-family:Georgia,'Times New Roman',serif;font-size:14px;color:var(--ink);white-space:nowrap;overflow:hidden;text-overflow:ellipsis;cursor:help}
.section{margin-top:28px}.section-head{display:flex;align-items:flex-end;justify-content:space-between;border-bottom:1px solid var(--ink);padding-bottom:7px;margin-bottom:14px}
.section-title{font-size:11px;font-weight:750;letter-spacing:.14em;text-transform:uppercase;color:var(--ink)}.section-note{font-size:10px;color:var(--muted)}
.metric-table{width:100%;border-collapse:collapse;background:#fff;border-bottom:1px solid var(--line);font-variant-numeric:tabular-nums}
.metric-table th{padding:8px 12px;border-bottom:1px solid var(--line);color:var(--muted);font-size:9px;text-transform:uppercase;letter-spacing:.11em;text-align:right;font-weight:650}
.metric-table th:first-child{text-align:left}.metric-table td{padding:10px 12px;border-bottom:1px solid #eceff1;text-align:right;font-size:12px}
.metric-table td:first-child{text-align:left;font-weight:600}.metric-table tr:last-child td{border-bottom:0}
.figure-panel{background:#fff;border:1px solid var(--line);padding:12px 12px 4px}
.notes-box{background:#fff;border-left:3px solid var(--bronze);padding:11px 15px;font-size:11px;color:#3e4b54;line-height:1.55}.notes-box p{margin:3px 0}
.dataTables_wrapper{font-size:11px}table.dataTable{border-collapse:collapse!important}
table.dataTable thead th{background:#f4f5f5!important;color:#4e5961!important;border-bottom:1px solid #cbd1d5!important;font-size:9px!important;text-transform:uppercase;letter-spacing:.06em}
table.dataTable tbody td{border-bottom:1px solid #eceff1!important}.dataTables_info,.dataTables_paginate{font-size:10px!important;color:var(--muted)!important}
@media(max-width:900px){.research-shell{grid-template-columns:1fr}.research-sidebar{border-right:0;border-bottom:1px solid var(--line)}.trial-strip{grid-template-columns:repeat(2,1fr)}.workspace{padding:22px 16px 40px}}

# ---- TALOS restrained single-trial workspace ----
.trial-line{
  font-family:Georgia,'Times New Roman',serif;
  font-size:16px;
  color:var(--ink);
  margin:5px 0 30px;
  padding-bottom:13px;
  border-bottom:1px solid var(--line);
}
.trial-line .dot{color:#9a7444;padding:0 8px;}
.research-section{margin:0 0 34px;}
.research-section-title{
  font-size:11px;font-weight:750;letter-spacing:.15em;text-transform:uppercase;
  color:var(--ink);border-bottom:1px solid var(--ink);padding-bottom:7px;margin-bottom:14px;
}
.compact-actions{
  margin-top:30px;padding-top:13px;border-top:1px solid var(--line);
  display:flex;align-items:center;gap:14px;flex-wrap:wrap;
}
.qc-pass{
  font-size:11px;font-weight:700;color:#526d59;letter-spacing:.03em;
}
.text-action{
  background:transparent!important;border:0!important;box-shadow:none!important;
  color:var(--navy)!important;font-size:11px!important;padding:4px 0!important;
  text-decoration:underline;text-underline-offset:3px;
}
.text-action:hover{color:var(--bronze)!important;background:transparent!important;}
.collapse-panel{
  margin-top:12px;background:#fff;border:1px solid var(--line);padding:12px;
}
.footer-disclaimer{
  margin-top:34px;padding-top:12px;border-top:1px solid #e5e8ea;
  font-size:9px;color:#7b858c;
}
"



app_css <- paste0(app_css, "
body{background:#f7f7f5!important;}
.research-header{
  height:86px!important;
  padding:0 48px!important;
  box-shadow:none!important;
}
.wordmark{font-size:30px!important;letter-spacing:.16em!important;}
.wordmark-sub{font-size:9px!important;letter-spacing:.17em!important;}
.header-rule{width:108px!important;height:2px!important;margin-top:8px!important;}

.research-shell{
  grid-template-columns:230px minmax(0,1fr)!important;
}
.research-sidebar{
  padding:30px 22px!important;
  background:#f1f3f3!important;
}
.side-title{margin-bottom:24px!important;}
.upload-box{margin-bottom:24px!important;}
.control-label{
  font-size:10px!important;
  letter-spacing:.03em!important;
  margin-bottom:9px!important;
}
.input-group{width:100%!important;}
.input-group-btn .btn{
  height:34px!important;
  padding:7px 10px!important;
}
.form-control{
  height:34px!important;
  background:#fff!important;
}
.progress{
  height:2px!important;
  margin-top:5px!important;
  background:#dce1e3!important;
  border-radius:0!important;
  box-shadow:none!important;
}
.progress-bar{
  background:#9a7444!important;
  box-shadow:none!important;
}
.status{
  margin-top:26px!important;
  background:transparent!important;
  padding:7px 0 7px 12px!important;
}

.workspace{
  max-width:1120px!important;
  padding:38px 54px 70px!important;
}
.trial-line{
  font-size:15px!important;
  line-height:1.5!important;
  margin:0 0 34px!important;
  padding-bottom:16px!important;
  border-bottom:1px solid #cfd5d8!important;
}
.trial-line .dot{padding:0 10px!important;color:#9a7444!important;}

.research-section{
  margin-bottom:44px!important;
}
.research-section-title{
  font-size:10px!important;
  letter-spacing:.18em!important;
  padding-bottom:9px!important;
  margin-bottom:18px!important;
}

.metric-table{
  max-width:820px!important;
  border:1px solid #d8dde1!important;
}
.metric-table th{
  background:#f5f6f6!important;
  padding:10px 16px!important;
}
.metric-table td{
  padding:13px 16px!important;
  font-size:12px!important;
}

.figure-panel{
  border:1px solid #d8dde1!important;
  padding:18px 18px 10px!important;
  box-shadow:none!important;
}
.js-plotly-plot .plotly .modebar{
  opacity:.18!important;
}
.js-plotly-plot:hover .plotly .modebar{
  opacity:.65!important;
}
.compact-actions{
  margin-top:38px!important;
  padding-top:16px!important;
}
.footer-disclaimer{margin-top:26px!important;}

@media(max-width:1100px){
  .workspace{padding:30px 28px 60px!important;}
}
")



safe_filename <- function(x) {
  x <- safe_chr(x)
  if (is.na(x) || !nzchar(x)) x <- "trial"
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "trial"
  x
}

talos_html_escape <- function(x) {
  htmltools::htmlEscape(as.character(x))
}

talos_timing_report_table <- function(r) {
  t <- r$timing_summary

  left_stride <- get_metric(t, "Stride time", "Left")
  right_stride <- get_metric(t, "Stride time", "Right")
  stride_asym <- get_metric(t, "Stride time", "Asymmetry")

  left_swing <- get_metric(t, "Swing time", "Left")
  right_swing <- get_metric(t, "Swing time", "Right")
  swing_asym <- get_metric(t, "Swing time", "Asymmetry")

  left_stance <- get_metric(t, "Stance time", "Left")
  right_stance <- get_metric(t, "Stance time", "Right")
  stance_asym <- get_metric(t, "Stance time", "Asymmetry")

  row_html <- function(label, left, right, asymmetry) {
    paste0(
      "<tr>",
      "<td>", talos_html_escape(label), "</td>",
      "<td>", talos_html_escape(fmt_value(left, 3, " s")), "</td>",
      "<td>", talos_html_escape(fmt_value(right, 3, " s")), "</td>",
      "<td>", talos_html_escape(fmt_value(asymmetry, 2, "%")), "</td>",
      "</tr>"
    )
  }

  paste0(
    "<table>",
    "<thead><tr>",
    "<th>Measure</th><th>Left</th><th>Right</th><th>Asymmetry</th>",
    "</tr></thead>",
    "<tbody>",
    row_html("Stride duration", left_stride, right_stride, stride_asym),
    row_html("Swing duration", left_swing, right_swing, swing_asym),
    row_html("Stance duration", left_stance, right_stance, stance_asym),
    "</tbody></table>"
  )
}

talos_qc_report_html <- function(r) {
  warnings <- r$qc_warnings

  if (length(warnings) == 0) {
    return("<p class='qc-pass'>QC: Passed</p>")
  }

  items <- paste0(
    "<li>",
    vapply(warnings, talos_html_escape, character(1)),
    "</li>",
    collapse = ""
  )

  paste0(
    "<p class='qc-review'>QC: Review</p>",
    "<ul>", items, "</ul>"
  )
}





app_css <- paste0(app_css, "
body{background:#eef1f2!important;color:#17232d!important;}
.research-header{display:none!important;}
.clinical-page{
  max-width:1280px;
  margin:26px auto 60px;
  background:#fff;
  border:1px solid #cfd5d8;
  box-shadow:0 1px 5px rgba(22,40,52,.07);
  padding:30px 38px 38px;
}
.clinical-top{
  display:flex;align-items:flex-start;justify-content:space-between;
  padding-bottom:12px;border-bottom:1.5px solid #24343f;
}
.clinical-brand{
  font-family:Georgia,'Times New Roman',serif;
  font-size:20px;letter-spacing:.09em;color:#173047;
}
.clinical-brand span{font-family:Arial,sans-serif;font-size:11px;letter-spacing:.14em;color:#68747d;margin-left:9px;}
.clinical-date{font-size:9px;letter-spacing:.12em;text-transform:uppercase;color:#6e7880;padding-top:5px;}
.clinical-meta{
  display:flex;gap:0;flex-wrap:wrap;
  padding:10px 0 11px;border-bottom:1px solid #d7dcdf;
  font-size:10px;color:#43515b;
}
.meta-item{padding:0 13px;border-right:1px solid #d7dcdf;}
.meta-item:first-child{padding-left:0}.meta-item:last-child{border-right:0}
.meta-key{font-weight:700;color:#253640;margin-right:5px;}
.report-section{margin-top:22px;}
.report-heading{
  font-family:Georgia,'Times New Roman',serif;
  font-size:15px;font-weight:700;color:#263944;
  margin:0 0 9px;padding-bottom:5px;border-bottom:1px solid #7f8a91;
}
.report-grid{display:grid;grid-template-columns:1fr 1fr;gap:20px;}
.report-panel{min-width:0;}
.clinical-table{width:100%;border-collapse:collapse;font-size:10px;font-variant-numeric:tabular-nums;}
.clinical-table th{
  padding:6px 9px;background:#f4f6f6;border-top:1px solid #d7dcdf;border-bottom:1px solid #bfc7cb;
  text-align:right;font-size:8px;letter-spacing:.07em;text-transform:uppercase;color:#52616a;
}
.clinical-table th:first-child{text-align:left;}
.clinical-table td{padding:7px 9px;border-bottom:1px solid #e3e6e8;text-align:right;}
.clinical-table td:first-child{text-align:left;font-weight:600;}
.figure-clean{border:0;background:#fff;padding:2px 0 0;}
.report-footer{
  display:flex;justify-content:flex-end;align-items:center;gap:18px;
  margin-top:24px;padding-top:10px;border-top:1px solid #cfd5d8;
}
.report-link{
  background:transparent!important;border:0!important;box-shadow:none!important;
  padding:2px 0!important;color:#173047!important;font-family:Georgia,'Times New Roman',serif!important;
  font-size:10px!important;text-decoration:underline;text-underline-offset:3px;
}
.report-link:hover{color:#9a7444!important;background:transparent!important;}
.qc-inline{margin-right:auto;font-size:9px;color:#526d59;font-weight:700;letter-spacing:.05em;}
.load-trial-wrap{position:fixed;right:22px;bottom:20px;z-index:20;}
.load-trial-btn{
  border:1px solid #173047!important;background:#173047!important;color:#fff!important;
  border-radius:2px!important;font-size:10px!important;padding:7px 12px!important;
  box-shadow:0 1px 3px rgba(0,0,0,.12)!important;
}
.modal-content{border-radius:2px!important;border:1px solid #bcc5ca!important;box-shadow:0 8px 28px rgba(0,0,0,.16)!important;}
.modal-header{border-bottom:1px solid #d7dcdf!important;}
.modal-title{font-family:Georgia,'Times New Roman',serif!important;color:#173047!important;}
.status{font-size:10px!important;}
.cycle-drawer{margin-top:12px;border:1px solid #d7dcdf;padding:12px;background:#fff;}
.analysis-panel{margin-top:14px;border:1px solid #d7dcdf;background:#fff;}
.analysis-panel summary{list-style:none;cursor:pointer;padding:11px 14px;font-family:Georgia,'Times New Roman',serif;font-size:13px;font-weight:700;color:#263944;background:#f8f9f9;border-bottom:0;display:flex;align-items:center;justify-content:space-between;}
.analysis-panel summary::-webkit-details-marker{display:none;}
.analysis-panel summary:after{content:'+';font-family:Arial,sans-serif;font-size:16px;font-weight:400;color:#69757d;}
.analysis-panel[open] summary{border-bottom:1px solid #d7dcdf;}
.analysis-panel[open] summary:after{content:'−';}
.panel-body{padding:12px 14px 14px;}
.cycle-panel-body{padding:0!important;width:100%;overflow:hidden;}
.cycle-panel-body .dataTables_wrapper{width:100%!important;margin:0!important;padding:12px 14px 10px!important;}
.cycle-panel-body table.dataTable{width:100%!important;margin:0!important;}
.cycle-panel-body .dataTables_scroll,.cycle-panel-body .dataTables_scrollHead,.cycle-panel-body .dataTables_scrollBody{width:100%!important;}
.panel-note{font-size:9px;color:#6f7a81;margin:7px 0 0;line-height:1.45;}
.report-disclaimer{font-size:8px;color:#7b858c;margin-top:18px;}
.js-plotly-plot .plotly .modebar{opacity:.10!important;}
.js-plotly-plot:hover .plotly .modebar{opacity:.55!important;}
@media(max-width:900px){
  .clinical-page{margin:0;border:0;padding:22px 16px;box-shadow:none}
  .report-grid{grid-template-columns:1fr}
  .clinical-meta{display:block}.meta-item{border-right:0;padding:3px 0}
}
")

ui <- fluidPage(
  tags$head(
    tags$title("TALOS / Gait Analysis"),
    tags$style(HTML(app_css))
  ),

  uiOutput("clinical_ui"),

  div(
    class = "load-trial-wrap",
    actionButton(
      "open_loader",
      "Load trial",
      class = "load-trial-btn"
    )
  )
)

server <- function(input, output, session) {

  result <- reactiveVal(NULL)
  analysis_error <- reactiveVal(NULL)
  analyzing <- reactiveVal(FALSE)

  observeEvent(input$open_loader, {
    showModal(
      modalDialog(
        title = "Load TALOS trial",

        fileInput(
          "processed_file",
          "Processed IMU data",
          multiple = TRUE,
          accept = ".txt",
          buttonLabel = "Select file(s)",
          placeholder = "No file selected"
        ),

        fileInput(
          "metadata_file",
          "Trial metadata",
          accept = c(".json", "application/json"),
          buttonLabel = "Select metadata",
          placeholder = "No file selected"
        ),

        uiOutput("loader_status"),

        easyClose = TRUE,
        footer = modalButton("Close")
      )
    )
  })

  observe({
    req(input$processed_file, input$metadata_file)

    analyzing(TRUE)
    analysis_error(NULL)

    tryCatch({
      r <- analyze_trial_hardened(
        input$processed_file$datapath,
        input$metadata_file$datapath,
        processed_names = input$processed_file$name
      )

      result(r)
      removeModal()

    }, error = function(e) {
      result(NULL)
      analysis_error(conditionMessage(e))
    })

    analyzing(FALSE)
  })

  output$loader_status <- renderUI({
    if (isTRUE(analyzing())) {
      return(div(class = "status wait", "Analyzing trial…"))
    }

    if (!is.null(analysis_error())) {
      return(
        div(
          class = "status error",
          paste("Could not analyze:", analysis_error())
        )
      )
    }

    NULL
  })

  variability_row <- function(label, left, right) {
    fmt <- function(x) if (length(x) == 0 || !is.finite(x)) "—" else sprintf("%.2f%%", x)
    tags$tr(tags$td(label), tags$td(fmt(left)), tags$td(fmt(right)))
  }

  count_row <- function(label, left, right) {
    fmt <- function(x) if (length(x) == 0 || !is.finite(x)) "—" else as.character(as.integer(x))
    tags$tr(tags$td(label), tags$td(fmt(left)), tags$td(fmt(right)))
  }

  output$clinical_ui <- renderUI({

    r <- result()

    if (is.null(r)) {
      return(
        div(
          class = "clinical-page",
          div(
            class = "clinical-top",
            div(class = "clinical-brand", "TALOS", tags$span("/ GAIT ANALYSIS")),
            div(class = "clinical-date", format(Sys.Date(), "%d %B %Y"))
          ),
          div(
            style = "padding:90px 20px;text-align:center;color:#71808a;",
            tags$div(
              style = "font-family:Georgia,serif;font-size:17px;color:#304550;margin-bottom:8px;",
              "No trial loaded"
            ),
            tags$div(
              style = "font-size:10px;",
              "Use Load trial to select processed IMU data and matching metadata."
            )
          )
        )
      )
    }

    m <- r$metadata
    t <- r$timing_summary

    val <- function(x, default = "—") {
      z <- safe_chr(x)
      if (is.na(z) || z == "") default else z
    }

    subject <- val(m$subject, "Uploaded trial")
    pathology <- val(m$pathology)
    deficit <- val(m$clinicalDeficitSide)

    left_stride <- get_metric(t, "Stride time", "Left")
    right_stride <- get_metric(t, "Stride time", "Right")
    stride_asym <- get_metric(t, "Stride time", "Asymmetry")

    left_swing <- get_metric(t, "Swing time", "Left")
    right_swing <- get_metric(t, "Swing time", "Right")
    swing_asym <- get_metric(t, "Swing time", "Asymmetry")

    left_stance <- get_metric(t, "Stance time", "Left")
    right_stance <- get_metric(t, "Stance time", "Right")
    stance_asym <- get_metric(t, "Stance time", "Asymmetry")

    timing_row <- function(label, left, right, asymmetry) {
      tags$tr(
        tags$td(label),
        tags$td(fmt_value(left, 3, " s")),
        tags$td(fmt_value(right, 3, " s")),
        tags$td(fmt_value(asymmetry, 2, "%"))
      )
    }

    tagList(
      div(
        class = "clinical-page",

        div(
          class = "clinical-top",
          div(
            class = "clinical-brand",
            "TALOS",
            tags$span("/ GAIT ANALYSIS")
          ),
          div(
            class = "clinical-date",
            toupper(format(Sys.Date(), "%d %B %Y"))
          )
        ),

        div(
          class = "clinical-meta",

          div(
            class = "meta-item",
            tags$span(class = "meta-key", "Subject:"),
            subject
          ),

          div(
            class = "meta-item",
            tags$span(class = "meta-key", "Condition:"),
            tags$span(title = pathology, pathology)
          ),

          div(
            class = "meta-item",
            tags$span(class = "meta-key", "Clinical deficit:"),
            tools::toTitleCase(deficit)
          ),

          div(
            class = "meta-item",
            tags$span(class = "meta-key", "QC:"),
            if (length(r$qc_warnings) == 0) "Passed" else "Review"
          )
        ),

        tags$details(
          class = "analysis-panel",
          open = NA,
          tags$summary("1. Temporal Gait Analysis"),
          div(
            class = "panel-body",
            tags$table(
              class = "clinical-table",
              tags$thead(tags$tr(tags$th("Measure"), tags$th("Left"), tags$th("Right"), tags$th("Asymmetry"))),
              tags$tbody(
                timing_row("Stride duration", left_stride, right_stride, stride_asym),
                timing_row("Swing duration", left_swing, right_swing, swing_asym),
                timing_row("Stance duration", left_stance, right_stance, stance_asym)
              )
            )
          )
        ),

        tags$details(
          class = "analysis-panel",
          tags$summary("2. Gait Variability & Consistency"),
          div(
            class = "panel-body",
            tags$table(
              class = "clinical-table",
              tags$thead(tags$tr(tags$th("Measure"), tags$th("Left"), tags$th("Right"))),
              tags$tbody(
                variability_row("Stride variability (CV)", r$left_summary$stride_cv_percent, r$right_summary$stride_cv_percent),
                variability_row("Swing variability (CV)", r$left_summary$swing_cv_percent, r$right_summary$swing_cv_percent),
                variability_row("Stance variability (CV)", r$left_summary$stance_cv_percent, r$right_summary$stance_cv_percent),
                count_row("Analyzed complete strides", r$left_summary$n_complete_strides, r$right_summary$n_complete_strides)
              )
            ),
            div(
              class = "panel-note",
              "CV = standard deviation / mean × 100. Lower values indicate more consistent cycle timing within this trial; these are descriptive measures, not diagnostic thresholds."
            )
          )
        ),

        tags$details(
          class = "analysis-panel",
          tags$summary("3. Bilateral Gait Events"),
          div(class = "panel-body figure-clean", plotlyOutput("gait_plot", height = "430px"))
        ),

        tags$details(
          class = "analysis-panel",
          tags$summary("4. Segmental IMU Motion"),
          div(class = "panel-body figure-clean", plotlyOutput("sensor_plot", height = "430px"))
        ),

        tags$details(
          class = "analysis-panel",
          tags$summary("5. Cycle-Level Data"),
          div(class = "panel-body cycle-panel-body", DTOutput("cycle_table"))
        ),

        div(
          class = "report-footer",

          uiOutput("qc_inline"),

          downloadButton(
            "download_report",
            "Export Report",
            class = "report-link"
          )
        ),

        div(
          class = "report-disclaimer",
          "Research-use software. Descriptive gait metrics are not diagnostic thresholds."
        )
      )
    )
  })

  output$qc_inline <- renderUI({
    req(result())

    w <- result()$qc_warnings

    tags$span(
      class = "qc-inline",
      title = if (length(w) == 0) {
        "Input validation completed without QC warnings."
      } else {
        paste(w, collapse = " | ")
      },
      if (length(w) == 0) "QC: ● Passed" else paste0("QC: Review (", length(w), ")")
    )
  })

  output$gait_plot <- renderPlotly({
    req(result())

    p <- make_gait_plot(result())

    plotly::config(
      p,
      displaylogo = FALSE,
      responsive = TRUE,
      modeBarButtonsToRemove = c(
        "select2d",
        "lasso2d",
        "autoScale2d",
        "toggleSpikelines",
        "hoverClosestCartesian",
        "hoverCompareCartesian"
      )
    )
  })

  output$sensor_plot <- renderPlotly({
    req(result())

    p <- make_sensor_plot(result())

    plotly::config(
      p,
      displaylogo = FALSE,
      responsive = TRUE,
      modeBarButtonsToRemove = c(
        "select2d",
        "lasso2d",
        "autoScale2d",
        "toggleSpikelines",
        "hoverClosestCartesian",
        "hoverCompareCartesian"
      )
    )
  })

  output$cycle_table <- renderDT({
    req(result())

    tab <- result()$cycles

    numeric_cols <- intersect(
      c("swing_s", "stance_s", "stride_s"),
      names(tab)
    )

    tab[numeric_cols] <- lapply(
      tab[numeric_cols],
      function(x) round(x, 3)
    )

    display_names <- c(
      side = "Side",
      phase = "Walking phase",
      cycle = "Cycle",
      event_number = "Event",
      swing_start = "Swing start",
      swing_end = "Swing end",
      swing_s = "Swing (s)",
      stance_s = "Stance (s)",
      stride_s = "Stride (s)"
    )

    keep <- intersect(names(display_names), names(tab))
    tab <- tab[, keep, drop = FALSE]
    names(tab) <- unname(display_names[keep])

    datatable(
      tab,
      rownames = FALSE,
      filter = "none",
      class = "stripe hover",
      options = list(
        pageLength = 12,
        scrollX = TRUE,
        autoWidth = FALSE,
        width = "100%",
        dom = "tip",
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      )
    )
  })

  output$download_report <- downloadHandler(

    filename = function() {
      r <- result()
      req(r)

      paste0(
        "TALOS_",
        safe_filename(safe_chr(r$metadata$subject)),
        "_gait_report.html"
      )
    },

    contentType = "text/html",

    content = function(file) {

      r <- result()
      req(r)

      subject <- safe_chr(r$metadata$subject)
      pathology <- safe_chr(r$metadata$pathology)
      deficit <- safe_chr(r$metadata$clinicalDeficitSide)
      fs <- safe_num(r$metadata$freq)

      if (is.na(subject) || !nzchar(subject)) subject <- "Uploaded trial"
      if (is.na(pathology) || !nzchar(pathology)) pathology <- "Not specified"
      if (is.na(deficit) || !nzchar(deficit)) deficit <- "Not specified"

      report_date <- format(Sys.Date(), "%d %B %Y")

      qc_html <- talos_qc_report_html(r)
      timing_html <- talos_timing_report_table(r)

      file_count <- if (!is.null(r$processed_file_count)) r$processed_file_count else 1L
      sample_count <- if (!is.null(r$data)) nrow(r$data) else NA_integer_

      html <- paste0(
        "<!doctype html>",
        "<html><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>TALOS Gait Analysis — ", talos_html_escape(subject), "</title>",
        "<style>",
        "body{font-family:Arial,Helvetica,sans-serif;background:#f3f5f5;color:#18232d;margin:0;padding:32px;}",
        ".page{max-width:950px;margin:auto;background:white;border:1px solid #cfd5d8;padding:34px 40px;}",
        ".top{display:flex;justify-content:space-between;border-bottom:1.5px solid #263944;padding-bottom:12px;}",
        ".brand{font-family:Georgia,'Times New Roman',serif;font-size:24px;letter-spacing:.08em;color:#173047;}",
        ".date{font-size:10px;letter-spacing:.08em;color:#6f7b83;text-transform:uppercase;padding-top:6px;}",
        ".meta{font-size:12px;padding:13px 0;border-bottom:1px solid #d7dcdf;line-height:1.8;}",
        "h2{font-family:Georgia,'Times New Roman',serif;font-size:16px;margin-top:26px;border-bottom:1px solid #7f8a91;padding-bottom:6px;color:#263944;}",
        "table{border-collapse:collapse;width:100%;font-size:12px;font-variant-numeric:tabular-nums;}",
        "th{background:#f4f6f6;color:#52616a;font-size:10px;text-transform:uppercase;letter-spacing:.06em;}",
        "th,td{padding:9px 10px;border-bottom:1px solid #e0e4e6;text-align:right;}",
        "th:first-child,td:first-child{text-align:left;}",
        ".qc-pass{color:#526d59;font-weight:700;}",
        ".qc-review{color:#956b2f;font-weight:700;}",
        ".methods{font-size:11px;line-height:1.65;color:#4d5a63;}",
        ".footer{margin-top:28px;padding-top:10px;border-top:1px solid #d7dcdf;font-size:9px;color:#7b858c;}",
        "</style></head><body><div class='page'>",

        "<div class='top'>",
        "<div class='brand'>TALOS / GAIT ANALYSIS</div>",
        "<div class='date'>", talos_html_escape(report_date), "</div>",
        "</div>",

        "<div class='meta'>",
        "<b>Subject:</b> ", talos_html_escape(subject),
        " &nbsp;&nbsp; | &nbsp;&nbsp; ",
        "<b>Condition:</b> ", talos_html_escape(pathology),
        " &nbsp;&nbsp; | &nbsp;&nbsp; ",
        "<b>Clinical deficit:</b> ", talos_html_escape(tools::toTitleCase(deficit)),
        "</div>",

        "<h2>1. Temporal Gait Analysis</h2>",
        timing_html,

        "<h2>2. Quality Control</h2>",
        qc_html,

        "<div class='methods'>",
        "<p><b>Sampling rate:</b> ", talos_html_escape(fmt_value(fs, 0, " Hz")), "<br>",
        "<b>Processed IMU file(s):</b> ", file_count, "<br>",
        "<b>Combined samples:</b> ", sample_count, "</p>",
        "<p>Cycle summaries are calculated from annotated straight-walking gait events; ",
        "cycles spanning the annotated U-turn are excluded from straight-walking timing summaries.</p>",
        "</div>",

        "<div class='footer'>",
        "Research-use software. Descriptive gait metrics are not diagnostic thresholds.",
        "</div>",

        "</div></body></html>"
      )

      con <- file(file, open = "wb")
      on.exit(close(con), add = TRUE)
      writeChar(enc2utf8(html), con, eos = NULL, useBytes = TRUE)
    }
  )
}

shinyApp(ui, server)
