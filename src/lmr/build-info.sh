#! /bin/bash

# Create json file of info
jq --null-input \
  --arg revisionId "$(git rev-parse --short HEAD)" \
  --arg fullRevisionId "$(git rev-parse HEAD)" \
  --arg revisionDate "$(git log -1 --format=%cd)" \
  --arg revisionAge "$(git log -1 --format=%cd --date=relative)" \
  --arg buildDate "$(date)" \
  '$ARGS.named' > build-info.json
