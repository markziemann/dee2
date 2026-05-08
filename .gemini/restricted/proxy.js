#!/usr/bin/env node

/**
 * Custom Proxy for Gemini Sandbox. Aim is to:
 * 1. Deny the ability to autheticate with GCP and make potentially destructive calls.
 * 2. Allow network traffic to minimally configured domains, by default those required
 *   to install codebase depenencies.
 *
 * Note: this is used partially because it can be easily installed and run within the
 * default Gemini sandbox. Other proxies could likely be adapted for a simliar job, but
 * may require extending the default Gemini sandbox docker image. This proxy can also
 * be fit to purpose.
 *
 * See https://geminicli.com/docs/examples/proxy-script/ for reference.
 */

import http from "node:http";
import net from "node:net";
import { URL } from "node:url";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const PROXY_PORT = 8877;
const ALLOWLIST_PATH = path.join(__dirname, "allowlist.txt");

// IPs/domains for the GCP Metadata server that we intercept and forward to our stub
const METADATA_HOSTS = ["169.254.169.254", "metadata.google.internal"];

let allowedDomains = [];
try {
  const content = fs.readFileSync(ALLOWLIST_PATH, "utf-8");
  allowedDomains = content
    .split("\n")
    .map((line) => line.trim())
    .filter((line) => line.length > 0 && !line.startsWith("#"));
} catch (err) {
  console.error(`[PROXY] Error reading allowlist.txt: ${err.message}`);
  allowedDomains = ["googleapis.com", "generativelanguage.googleapis.com"];
}

const server = http.createServer((req, res) => {
  // 1. Intercept GET requests for GCP Metadata server and forward to our stub server
  try {
    const url = new URL(req.url);
    if (METADATA_HOSTS.includes(url.hostname)) {
      console.log(`[PROXY] Intercepting GET (GCP Metadata) for: ${req.url}`);

      const STUB_PORT = process.env.STUB_METADATA_PORT || 8878;
      const STUB_HOST = "host.docker.internal";

      const proxyReq = http.request(
        {
          host: STUB_HOST,
          port: STUB_PORT,
          path: url.pathname + url.search,
          method: req.method,
          headers: req.headers,
        },
        (proxyRes) => {
          res.writeHead(proxyRes.statusCode, proxyRes.headers);
          proxyRes.pipe(res);
        },
      );

      proxyReq.on("error", (err) => {
        console.error(
          `[PROXY] Error forwarding GET to stub metadata server: ${err.message}`,
        );
        res.writeHead(502);
        res.end("Bad Gateway");
      });

      req.pipe(proxyReq);
      return;
    }
  } catch (err) {
    // Not a full URL in req.url, likely just a path. Ignore and fall through.
  }

  console.log(
    `[PROXY] Denying non-CONNECT request for: ${req.method} ${req.url}`,
  );
  res.writeHead(405, { "Content-Type": "text/plain" });
  res.end("Method Not Allowed");
});

server.on("connect", (req, clientSocket, head) => {
  let port, hostname;
  try {
    const url = new URL(`http://${req.url}`);
    port = url.port || "443";
    hostname = url.hostname;
  } catch (err) {
    console.error(`[PROXY] Invalid URL: ${req.url}`);
    clientSocket.end("HTTP/1.1 400 Bad Request\r\n\r\n");
    return;
  }

  // 1. Intercept GCP Metadata server and redirect to our stub server
  if (METADATA_HOSTS.includes(hostname)) {
    console.log(
      `[PROXY] REDIRECTING (GCP Metadata) to STUB: ${hostname}:${port}`,
    );

    // Redirect metadata requests to our stub server running on the host
    const STUB_PORT = process.env.STUB_METADATA_PORT || 8878;
    const STUB_HOST = "host.docker.internal";

    const serverSocket = net.connect(STUB_PORT, STUB_HOST, () => {
      clientSocket.write("HTTP/1.1 200 Connection Established\r\n\r\n");
      serverSocket.write(head);
      serverSocket.pipe(clientSocket);
      clientSocket.pipe(serverSocket);
    });

    serverSocket.on("error", (err) => {
      console.error(
        `[PROXY] Error connecting to stub metadata server: ${err.message}`,
      );
      clientSocket.end(`HTTP/1.1 502 Bad Gateway\r\n\r\n`);
    });

    return;
  }

  // 2. Check allowlist (loaded on startup)
  const isAllowed = allowedDomains.some(
    (domain) => hostname === domain || hostname.endsWith(`.${domain}`),
  );

  if (isAllowed) {
    // To debug all traffic, add: console.log(`[PROXY] ALLOWED: ${hostname}:${port}`);
    const serverSocket = net.connect(port, hostname, () => {
      clientSocket.write("HTTP/1.1 200 Connection Established\r\n\r\n");
      serverSocket.write(head);
      serverSocket.pipe(clientSocket);
      clientSocket.pipe(serverSocket);
    });

    serverSocket.on("error", (err) => {
      console.error(
        `[PROXY] Error connecting to destination ${hostname}: ${err.message}`,
      );
      clientSocket.end(`HTTP/1.1 502 Bad Gateway\r\n\r\n`);
    });
  } else {
    console.log(`[PROXY] DENIED (Add to allowlist.txt): ${hostname}`);
    clientSocket.end("HTTP/1.1 403 Forbidden\r\n\r\n");
  }

  clientSocket.on("error", (err) => {
    // Silently handle client socket hangs
  });
});

server.listen(PROXY_PORT, () => {
  const address = server.address();
  console.log(
    `[PROXY] Sandbox Proxy listening on ${address.address}:${address.port}`,
  );
  console.log(`[PROXY] Loaded static allowlist from ${ALLOWLIST_PATH}`);
});
