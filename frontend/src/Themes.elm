module Themes exposing (..)

import Html exposing (Attribute)
import Html.Attributes exposing (..)


type alias AttributeUpdater msg =
    List (Attribute msg) -> List (Attribute msg)


type alias Theme msg =
    { background : AttributeUpdater msg
    , textColour : AttributeUpdater msg
    }


defaultTheme =
    { background = addStyle (Just <| class "bg-light")
    , textColour = addStyle (Just <| class "text-primary")
    }


type alias NavTheme msg =
    Theme msg


defaultNavTheme : Theme msg
defaultNavTheme =
    { defaultTheme
        | background = addStyle (Just <| style "background-color" "rgb(204, 230, 245)")
            , textColour =
                addStyle
                    (Just <| style "color" "rgb(5, 86, 122)")
    }


type alias SiteTheme msg =
    { navTheme : NavTheme msg }


defaultSiteTheme =
    { navTheme = defaultNavTheme }


addStyle : Maybe (Attribute msg) -> List (Attribute msg) -> List (Attribute msg)
addStyle maybeStyle attributes =
    case maybeStyle of
        Just style ->
            style :: attributes

        Nothing ->
            attributes
