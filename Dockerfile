FROM rust:1.49.0 as builder
WORKDIR /usr/src/methylcx
COPY . .
RUN cargo install --path .

FROM debian:buster-slim
COPY --from=builder /usr/local/cargo/bin/methylcx /usr/local/bin/methylcx
CMD ["methylcx"]