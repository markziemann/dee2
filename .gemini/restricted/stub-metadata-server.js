/**
 * Stub GCE Metadata Server for Gemini Restricted Sandbox
 *
 * This server runs on the host machine and simulates the Google Compute Engine
 * (GCE) metadata service (normally accessible at 169.254.169.254).
 *
 * Why is this needed?
 * We want the Gemini CLI sandbox to authenticate to GCP services using a highly
 * restricted impersonated service account, rather than the host's personal
 * credentials (Application Default Credentials).
 *
 * To accomplish this without exposing the host's raw gcloud binaries into the sandbox:
 * 1. The sandbox's network proxy intercepts requests to `metadata.google.internal`
 *    and redirects them to this local Node.js server.
 * 2. When the Google Auth SDK inside the sandbox fails to find local credential
 *    files, it falls back to querying the metadata server for tokens.
 * 3. This script fields those token requests by executing `gcloud auth print-access-token
 *    --impersonate-service-account` on the host, caching the result, and serving it
 *    back to the sandbox in the format the GCP SDKs expect.
 */
import http from "node:http";
import { exec } from "node:child_process";
import util from "node:util";

const execAsync = util.promisify(exec);

const PORT = process.env.STUB_METADATA_PORT || 8878;
const IMPERSONATE_SA =
  process.env.IMPERSONATE_SA ||
  "gemini-cli-restricted@x-alife-analysis.iam.gserviceaccount.com";
const GOOGLE_PROJECT = process.env.GOOGLE_PROJECT || "x-alife-analysis";

let cachedToken = null;
let tokenExpiry = 0; // Epoch ms

async function getAccessToken() {
  const now = Date.now();
  if (cachedToken && now < tokenExpiry) {
    return cachedToken;
  }

  console.log(
    `[STUB METADATA] Fetching new token via gcloud for ${IMPERSONATE_SA}...`,
  );
  try {
    const { stdout } = await execAsync(
      `gcloud auth print-access-token --impersonate-service-account=${IMPERSONATE_SA}`,
    );
    const token = stdout.trim();
    cachedToken = token;
    // Cache for 45 minutes (tokens usually live for 60m)
    tokenExpiry = now + 45 * 60 * 1000;
    return token;
  } catch (error) {
    console.error("[STUB METADATA] Error getting token:", error);
    throw error;
  }
}

const server = http.createServer(async (req, res) => {
  const timestamp = new Date().toISOString();
  console.log(
    `[${timestamp}] [STUB METADATA] Request: ${req.method} ${req.url}`,
  );
  console.log(
    `[${timestamp}] [STUB METADATA] Headers: ${JSON.stringify(req.headers)}`,
  );

  // Require the flavor header, per standard GCP metadata server
  if (req.headers["metadata-flavor"] !== "Google") {
    console.warn(
      `[${timestamp}] [STUB METADATA] Missing Metadata-Flavor header from ${req.socket.remoteAddress}`,
    );
    res.writeHead(403);
    res.end("Missing Metadata-Flavor: Google header");
    return;
  }

  try {
    if (req.url === "/computeMetadata/v1/project/project-id") {
      console.log(
        `[${timestamp}] [STUB METADATA] Serving project-id: ${GOOGLE_PROJECT}`,
      );
      res.writeHead(200, {
        "Content-Type": "application/text",
        "Metadata-Flavor": "Google",
      });
      res.end(GOOGLE_PROJECT);
    } else if (
      req.url.startsWith(
        "/computeMetadata/v1/instance/service-accounts/default/token",
      )
    ) {
      console.log(`[${timestamp}] [STUB METADATA] Serving access token...`);
      const token = await getAccessToken();
      const response = {
        access_token: token,
        expires_in: Math.floor((tokenExpiry - Date.now()) / 1000),
        token_type: "Bearer",
      };
      res.writeHead(200, {
        "Content-Type": "application/json",
        "Metadata-Flavor": "Google",
      });
      res.end(JSON.stringify(response));
    } else if (
      req.url === "/computeMetadata/v1/instance/service-accounts/default/email"
    ) {
      console.log(
        `[${timestamp}] [STUB METADATA] Serving email: ${IMPERSONATE_SA}`,
      );
      res.writeHead(200, {
        "Content-Type": "application/text",
        "Metadata-Flavor": "Google",
      });
      res.end(IMPERSONATE_SA);
    } else {
      console.log(`[${timestamp}] [STUB METADATA] 404 Not Found: ${req.url}`);
      res.writeHead(404, { "Metadata-Flavor": "Google" });
      res.end("Not Found");
    }
  } catch (err) {
    console.error(
      `[${timestamp}] [STUB METADATA] 500 Internal Error: ${err.message}`,
    );
    res.writeHead(500, { "Metadata-Flavor": "Google" });
    res.end("Internal Server Error");
  }
});

server.listen(PORT, "0.0.0.0", () => {
  console.log(`[STUB METADATA] Server listening on port ${PORT}`);
  console.log(`[STUB METADATA] Impersonating: ${IMPERSONATE_SA}`);
});
