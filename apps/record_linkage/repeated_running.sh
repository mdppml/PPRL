#!/bin/bash
#make sure the background processes are killed if the script is interrupted:
trap "kill 0" EXIT
# parse arguments:
JITTER=0
MODE=d
while [[ $# -gt 0 ]]; do
    i=$1
    case $i in
        -h|--help)
            help
            exit 0
            ;;
        -l=*|--latency=*)
            LATENCY="${i#*=}"
            shift
            ;;
        -j=*|--jitter=*)
            JITTER="${i#*=}"
            shift
            ;;
        -b=*|--bandwidth=*)
            BANDWIDTH="${i#*=}"
            shift
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1") # save positional arg
            shift
            ;;
    esac
done
# make sure that enough positional arguments were provided:
if [[ ${#POSITIONAL_ARGS[@]} -lt 2 ]]
then
  echo "Positional arguments were not provided. Call this tool with -h to see which ones are required."
  exit 1
fi
# make positional arguments into string:
ARGS=${POSITIONAL_ARGS[*]}
EXECUTABLE=../../../build/Release/record_linkage
# make sure that latency is larger than jitter:
if [[ "$JITTER" -gt 0 ]]
then
  if ! [[ "$LATENCY" ]] || [[ "$LATENCY" -lt "$JITTER" ]]
  then
    echo "Jitter cannot be higher than Latency."
    exit 1
  fi
fi
PORT_P1P2=6381
PORT_P2P1=6381
PORT_HELPERP=6379
PORT_PHELPER=6379
# set up toxiproxy if it is used:
if [[ "$LATENCY" ]] || [[ "$BANDWIDTH" ]]
then
  PORT_P2P1=26381
  PORT_PHELPER=26379
  # find toxiproxy executables:
  TOXI_SERVER_ARRAY=(toxiproxy-server*)
  TOXI_SERVER=${TOXI_SERVER_ARRAY[0]}
  TOXI_CLI_ARRAY=(toxiproxy-cli*)
  TOXI_CLI=${TOXI_CLI_ARRAY[0]}
  # initialise toxiproxy:
  "./$TOXI_SERVER" &>toxi.log &
  sleep 1
  "./$TOXI_CLI" create \
  -l localhost:"$((PORT_PHELPER))" \
  -u localhost:"$((PORT_HELPERP))" \
  "helper_p1" &>toxi.log
  "./$TOXI_CLI" create \
  -l localhost:"$((PORT_PHELPER + 1))" \
  -u localhost:"$((PORT_HELPERP + 1))" \
  "helper_p2" &>toxi.log
  "./$TOXI_CLI" create \
  -l localhost:"$((PORT_P2P1))" \
  -u localhost:"$((PORT_P1P2))" \
  "p1_p2" &>toxi.log

  if [[ "$LATENCY" ]]
  then
    echo "Added Latency: ${LATENCY}"
    echo "Jitter: ${JITTER}"
      "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --downstream "helper_p1" &>toxi.log
      "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --upstream "helper_p1" &>toxi.log
      "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --downstream "helper_p2" &>toxi.log
      "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --upstream "helper_p2" &>toxi.log
      "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --downstream "p1_p2" &>toxi.log
      "./$TOXI_CLI" toxic add -t latency -a latency="$LATENCY" -a jitter="$JITTER" --upstream "p1_p2" &>toxi.log
  fi
  if [[ "$BANDWIDTH" ]]
  then
    echo "Added Bandwidth: ${BANDWIDTH}"
      "./$TOXI_CLI" toxic add -t bandwidth -a rate="$BANDWIDTH" --downstream "helper_p1" &>toxi.log
      "./$TOXI_CLI" toxic add -t bandwidth -a rate="$BANDWIDTH" --upstream "helper_p1" &>toxi.log
      "./$TOXI_CLI" toxic add -t bandwidth -a rate="$BANDWIDTH" --downstream "helper_p2" &>toxi.log
      "./$TOXI_CLI" toxic add -t bandwidth -a rate="$BANDWIDTH" --upstream "helper_p2" &>toxi.log
      "./$TOXI_CLI" toxic add -t bandwidth -a rate="$BANDWIDTH" --downstream "p1_p2" &>toxi.log
      "./$TOXI_CLI" toxic add -t bandwidth -a rate="$BANDWIDTH" --upstream "p1_p2" &>toxi.log
  fi
fi
# run benchmarks:
for RUN in {0..10}
  do
  $EXECUTABLE helper $PORT_HELPERP 127.0.0.1 0 0 $ARGS &
  sleep 2
  $EXECUTABLE proxy1 $PORT_PHELPER 127.0.0.1 $PORT_P1P2 127.0.0.1 $ARGS &
  sleep 2
  $EXECUTABLE proxy2 $PORT_PHELPER 127.0.0.1 $PORT_P2P1 127.0.0.1 $ARGS > result${RUN}.csv
  done
