module ResultsPage.Views exposing (..)

import Array
import Bool.Extra as BE
import Dict
import Helpers exposing (errorToString)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events as Events exposing (onClick)
import ResultsPage.Helpers exposing (..)
import ResultsPage.Types exposing (..)
import SearchPage.Helpers exposing (withPagination)
import SearchPage.Types exposing (SearchParameters, SearchResult, SearchResults)
import Set
import SharedTypes exposing (PaginationOffset, RemoteData(..), WebData, unwrapWebData)
import Table exposing (HtmlDetails, Status)


selectClickedResult : (SearchResult -> Msg) -> SelectedResults -> SearchResult -> List (Html.Attribute Msg)
selectClickedResult msg selectedResults searchResult =
    [ class (highlightRowIfTrue "table-primary" <| Dict.member searchResult.id selectedResults)
    , style "cursor" "pointer" -- First place not using bootstrap for style?
    , Events.onClick (msg searchResult)
    ]


stageResultForRemoval : (Int -> Msg) -> ResultsPendingRemoval -> SelectedResult -> List (Html.Attribute Msg)
stageResultForRemoval msg resultsPendingRemoval ( id, ( _, _ ) ) =
    [ class (highlightRowIfTrue "table-danger" <| Set.member id resultsPendingRemoval)
    , style "cursor" "pointer" -- First place not using bootstrap for style?
    , Events.onClick (msg id)
    ]


resultsTable : SelectedResults -> Table.Config SearchResult Msg
resultsTable selectedResults =
    Table.customConfig
        { toId = getId
        , toMsg = SetResultsTableState
        , columns =
            [ Table.stringColumn "Row" getId
            , Table.stringColumn "Species" (get "species")
            , Table.stringColumn "SRA Run" (get "SRR_accession")
            , noOverflowColumn "QC summary" (\data -> ( data.id, get "QC_summary" data ))
            , Table.stringColumn "SRA experiment" (get "SRX_accession")
            , Table.stringColumn "SRA sample" (get "SRS_accession")
            , Table.stringColumn "SRA project" (get "SRP_accession")
            , noOverflowColumn "Sample" (\data -> ( data.id, get "Sample_name" data ))
            , Table.stringColumn "Experiment" (get "GEO_series")
            ]
        , customizations =
            { defaultTable
                | tableAttrs =
                    [ class "table table-hover table-sm table-bordered table-responsive-lg"
                    , style "table-layout" "fixed" -- Prevents table going wider than parent element
                    , style "font-size" "clamp(12px, 4vw, 14px)"
                    ]
                , rowAttrs = selectClickedResult ResultClicked selectedResults
            }
        }


selectedTable : ResultsPendingRemoval -> Table.Config SelectedResult Msg
selectedTable resultsPendingRemoval =
    Table.customConfig
        { toId = Tuple.first >> String.fromInt
        , toMsg = SetSelectedResultsTableState
        , columns =
            [ Table.intColumn "Row" Tuple.first
            , Table.stringColumn "Species" (Tuple.second >> Tuple.first)
            , Table.stringColumn "SRA" (Tuple.second >> Tuple.second)
            ]
        , customizations =
            { defaultTable
                | tableAttrs =
                    [ class "table table-hover table-sm table-bordered"
                    , style "table-layout" "fixed" -- Prevents table going wider than parent element
                    , style "font-size" "clamp(12px, 4vw, 14px)"
                    ]
                , rowAttrs = stageResultForRemoval SelectedResultClicked resultsPendingRemoval
            }
        }


buttonOrSpinner : Bool -> SelectedResults -> Html Msg
buttonOrSpinner downloading rows =
    if not downloading then
        a
            [ hideWhenTrue "btn btn-outline-primary btn-block" (Dict.isEmpty rows) |> class
            , href <| "api/download/" ++ queryString rows
            , attribute "download" "data.zip"
            , onClick DownloadRequested
            ]
            [ text "download" ]

    else
        button
            [ class "btn btn-primary btn-block"
            , attribute "type" "button"
            , attribute "disabled" "disabled"
            ]
            [ span
                [ class "spinner-border spinner-border-sm mr-2"
                , attribute "role" "status"
                , attribute "aria-hidden" "true"
                ]
                []
            , text "Loading..."
            ]


pageSelector disable page label =
    li [ class (disableWhenTrue "page-item" disable) ]
        [ button [ onClick (PageRequest page), class "page-link", attribute "tabindex" (BE.ifElse "-1" "0" disable) ]
            [ text label ]
        ]


pagination : SearchParameters -> PaginationOffset -> Int -> Html Msg
pagination searchParameters ({ perPage, offset } as paginationOffset) hits =
    let
        previousButton =
            if offset < perPage then
                -- Disabled state
                pageSelector
                    True
                    (withPagination paginationOffset searchParameters)
                    "Previous"

            else
                pageSelector
                    False
                    (withPagination (PaginationOffset perPage <| offset - perPage) searchParameters)
                    "Previous"

        nextButton =
            --Disabled state
            if (hits - offset) <= perPage then
                pageSelector
                    True
                    (withPagination paginationOffset searchParameters)
                    "Next"

            else
                pageSelector
                    False
                    (withPagination (PaginationOffset perPage <| offset + perPage) searchParameters)
                    "Next"
    in
    nav [ attribute "aria-label" "Page navigation" ]
        [ ul [ class "pagination justify-content-end" ]
            [ previousButton
            , nextButton
            ]
        ]


viewSearchResultHits : WebData SearchResults -> Html Msg
viewSearchResultHits searchResults =
    case searchResults of
        Success results ->
            text <| "Hits: " ++ String.fromInt results.hits

        Failure err ->
            div [ class "text-danger" ] [ text <| errorToString err ]

        _ ->
            div [ class "spinner-border spinner-border-sm", attribute "role" "status" ]
                [ span [ class "sr-only" ]
                    [ text "Loading..." ]
                ]


viewSearchResults : Model -> PaginationOffset -> List (Html Msg)
viewSearchResults ({ resultsTableState, resultsTableQuery } as model) paginationOffset =
    let
        maybeExpiredSearchResults =
            maybeExpiredData model.searchResults
    in
    case model.searchParameters of
        Just searchParameters ->
            [ div [ class "row" ]
                [ div [ class "col-xl-9" ]
                    [ div [ class "bg-light text-primary" ]
                        [ viewSearchResultHits maybeExpiredSearchResults ]
                    , Table.view
                        (resultsTable model.selectedResults)
                        resultsTableState
                        (unwrapWebData [] (.rows >> Array.toList) maybeExpiredSearchResults)
                    , pagination searchParameters paginationOffset (unwrapWebData 0 .hits maybeExpiredSearchResults)
                    ]
                , div [ class "col-xl-3" ]
                    [ div [ class "sticky-top" ]
                        [ div [ class "bg-light text-primary" ]
                            [ text <| "Selected: " ++ (Dict.size model.selectedResults |> String.fromInt) ]
                        , Table.view
                            (selectedTable model.resultsPendingRemoval)
                            model.selectedResultsTableState
                            (Dict.toList model.selectedResults)
                        , div []
                            [ buttonOrSpinner model.downloading model.selectedResults
                            , button
                                [ hideWhenTrue
                                    "btn btn-outline-danger btn-block"
                                    (Set.isEmpty model.resultsPendingRemoval)
                                    |> class
                                , onClick RemoveStagedSelections
                                ]
                                [ text "Remove" ]
                            ]
                        ]
                    ]
                ]
            ]

        Nothing ->
            [ div [] [ text "I need something to search for first!" ] ]
