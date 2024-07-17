# Reference Builder

Build and maintain reference sets of pathogen genome sequences.

## Installation

```shell script
pip install ref-builder
```

## Environmental Variables

Ref Builder uses the NCBI API to fetch genome data.

Unauthenticated requests are limited to 3 per second. Setting NCBI account credentials
in your environmental variables can increase the allowed rate to 10 requests per
second.

| Name | Description |
|----|---------|
| `NCBI_EMAIL` | The e-mail address used for your NCBI account |
| `NCBI_API_KEY` | The [API key](https://www.ncbi.nlm.nih.gov/account/settings/) associated with your NCBI account. |